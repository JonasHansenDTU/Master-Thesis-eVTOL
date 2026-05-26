function Destructor(plane::planeSolution, idxs::Vector{Int64}, data)
    m = Int(plane.flightLegs)

    if m == 2 && !isempty(idxs)
        base = plane.route[1]
        endpoint = plane.route[2]
        
        # Check if endpoint is valid single-leg endpoint (in end_vp and not base)
        end_VP_choices = data.end_vp[Int(base)]
        if endpoint in end_VP_choices && endpoint != base
            # Preserve as single-leg route
            resize!(plane.route, 2)
            empty!(plane.turnaroundTime)
            insert!(plane.turnaroundTime, 1, Int32(data.te))
            plane.flightLegs = Int32(1)
            return
        end
        
        # Otherwise reset to base (0 legs)
        resize!(plane.route, 1)
        plane.route[1] = base
        empty!(plane.turnaroundTime)
        plane.flightLegs = Int32(0)
        return
    end

    del = sort(unique(Int.(idxs)))

    # keep only internal stop indices for safety
    del = [idx for idx in del if 2 <= idx <= m]
    isempty(del) && return

    plane.flightLegs -= Int32(length(del))
    deleteat!(plane.route, del)
    deleteat!(plane.turnaroundTime, del)
end

function has_consecutive_duplicates(route)
    for i in 1:(length(route)-1)
        if route[i] == route[i+1]
            return true
        end
    end
    return false
end

function Constructor(plane::planeSolution, VP::Int, idx::Int, data)
    te32 = Int32(data.te)

    if plane.flightLegs == 0
        base = Int32(plane.route[1])
        end_VP_choices = data.end_vp[Int(base)]
        
        # Try to create single-leg route if VP is in allowed endpoints and not base
        if VP != base && VP in end_VP_choices
            insert!(plane.route, 2, Int32(VP))
            insert!(plane.turnaroundTime, 1, te32)
            plane.flightLegs = 1
            return
        end
        
        # Otherwise fall back to two-leg base-return route
        insert!(plane.route, 2, Int32(VP))
        insert!(plane.route, 3, base)

        insert!(plane.turnaroundTime, 1, te32)
        insert!(plane.turnaroundTime, 2, te32)

        plane.flightLegs += 2
        return
    end

    plane.flightLegs += 1
    insert!(plane.route, idx, Int32(VP))
    insert!(plane.turnaroundTime, idx, te32)

end

function two_opt(plane::planeSolution, idx1::Int, idx2::Int)

    # Nothing to improve if there are fewer than 2 flown legs.
    if plane.flightLegs < 2
        return plane
    end

    # Ensure ordered bounds.
    if idx1 > idx2
        idx1, idx2 = idx2, idx1
    end

    # Keep first/last node fixed (base), reverse only internal part.
    lo = max(2, idx1)
    hi = min(length(plane.route) - 1, idx2)

    if lo >= hi
        return plane
    end

    rev_seg = reverse(plane.route[lo:hi])

    # Reject move if it creates consecutive identical vertiports.
    if lo > 1 && plane.route[lo - 1] == rev_seg[1]
        return plane
    end
    if hi < length(plane.route) && rev_seg[end] == plane.route[hi + 1]
        return plane
    end
    for k in 1:(length(rev_seg) - 1)
        if rev_seg[k] == rev_seg[k + 1]
            return plane
        end
    end

    plane.route[lo:hi] = rev_seg
    return plane
end

function UpdateTurnAroundTimes(planes::allPlaneSolution, from::Int64, maxTurnaround::Int64, data)
    op = data.op
    dp = data.dp
    te = Int(data.te)
    dt = data.dt
    rt = data.rt
    w = data.w
    ET = data.ET

    for plane in planes.planes
        first_trip = (from == 1 ? true : false)
        flightLegs = Int(plane.flightLegs)
        if flightLegs == 0
            continue
        end

        route = plane.route
        from_k = clamp(Int(from), 1, flightLegs)

        # Rebuild state up to the first leg we want to update
        current_time = 0
        current_VP = Int(route[1])

        for k in 1:(from_k - 1)
            next_VP = Int(route[k + 1])
            current_time += Int(plane.turnaroundTime[k]) + rt[(current_VP, next_VP)]
            current_VP = next_VP
        end

        # Recompute turnaround times from from_k onward
        for k in from_k:flightLegs
            common_a = [a for (a, v) in op if v == current_VP && dp[a] == route[k+1]]

            candidate_pass = [dt[a] for a in common_a if
                current_time + (first_trip ? 0 : te) - w <= dt[a] <= current_time + maxTurnaround &&
                dt[a] + rt[(op[a], dp[a])] <= ET]

            if !isempty(candidate_pass) && rand() <= 0.5
                # deterministic: latest feasible passenger time in the window
                chosen_time = max(maximum(candidate_pass) - current_time, (first_trip ? 0 : te))
                plane.turnaroundTime[k] = Int32(round(Int, chosen_time))
            else
                # Choose rate (tune as needed)
                λ = 1.2

                # Bounds as floats
                a = float(first_trip ? 0 : te)
                b = float(maxTurnaround)

                # Sample from truncated exponential using inverse CDF
                u = rand()
                x = a - log(1 - u * (1 - exp(-λ * (b - a)))) / λ

                # Convert to integer safely
                x_int = floor(Int, x)

                plane.turnaroundTime[k] = x_int
            end

            next_VP = Int(route[k + 1])
            current_time += Int(plane.turnaroundTime[k]) + rt[(current_VP, next_VP)]
            current_VP = next_VP
            first_trip = false
        end
    end
end

function obj(planes::allPlaneSolution, data, rt)
    assignments, scheduled = assign_passengersV2(planes, data, Int.(rt))

    fitness = fitnessFunction(
            planes,
            assignments,
            Float32(data.bmax),
            Float32(data.bmid),
            Float32(data.bmin),
            data.dist,
            Float32(data.ec),
            Float32(data.battery_per_km),
            Int.(rt),
            Int(round(data.ET)),
            maximum(Int.(data.T)),
            maximum(data.V),
            data.cap_v,
            data
        )

    return fitness
end

function Filter(plane::planeSolution, rt, e, te, ec, ET, bmin, bmid, bmax)


    flightLegs = plane.flightLegs

    if flightLegs <= 1
        return true, false
    end

    Total_traveltime = 0
    Total_battery_needed = 0
    
    for i in 1:flightLegs
        Total_traveltime += rt[(plane.route[i], plane.route[i+1])]
        Total_battery_needed += e[(plane.route[i], plane.route[i+1])]
        if e[(plane.route[i], plane.route[i+1])] > (i == 1 ? bmid : bmax) - bmin        # Check if there are parts of the route using > bmax-bmin (or bmid-bmin if first trip)
            return false, true
        end
    end

    minimum_turnaround = (flightLegs-1)*te
    maximum_battery_charge = (ET-Total_traveltime)*ec


    Feasible_Turnaround = Total_traveltime + minimum_turnaround <= ET                   #Enough Time to complete the route?
    Feasible_Battery = bmid + maximum_battery_charge - Total_battery_needed >= bmin     #Feasible to charge battery in route?


     # TimeToCharge = ET - Total_traveltime >= (bmin-bmid + Total_battery_needed)*1/ec     #Connect the two and see if there can be time to charge and make it back (Turns out to be identical to Feasible_Battery)


    return (Feasible_Turnaround && Feasible_Battery) , false
end


function assign_Turnaround(Planes::allPlaneSolution, data)

    N = data.N
    bmin = data.bmin
    bmid = data.bmid
    bmax = data.bmax
    e = data.e
    ec = data.ec
    ET = data.ET
    te = data.te
    rt = data.rt

    for plane in Planes.planes

        #Only look at planes with infeasible routes


        first_trip = true
        m = plane.flightLegs
        CurrrentBattery = 80

        max_turnaorund = ET - (m-2)*te
                
        for i in 1:m
            max_turnaorund -= rt[(plane.route[i],plane.route[i+1])]
        end

        b_target = 80

        for i in 1:m

            if first_trip
                plane.turnaroundTime[i] = 0
                first_trip = false
            else
                
                max_turnaorund -= plane.turnaroundTime[i-1]
                
                CurrrentBattery -= e[(plane.route[i], plane.route[i+1])] 

                timeToBTarget = round((b_target - CurrrentBattery)/ec)

                plane.turnaroundTime[i] = max(min(max_turnaorund, timeToBTarget),te)

            end

        end

    end
    
end



function DestructLoop(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    N = data.N
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Do not delete first node (route[1]); last route node is not in this range anyway
        for idx in 2:m
            temp_sol = deepcopy(planes)

            Destructor(temp_sol.planes[n], [idx], data)

            if has_consecutive_duplicates(temp_sol.planes[n].route)
                Destructor(temp_sol.planes[n], [idx], data)
            end

            new_obj = obj(temp_sol, data, rt)

            temp_sol2 = deepcopy(temp_sol)
            from_idx = Int64(idx - 1)
            UpdateTurnAroundTimes(temp_sol2, from_idx, maxTurnaround, data)
            new_obj2 = obj(temp_sol2, data, rt)

            if new_obj2 > new_obj
                temp_sol = temp_sol2
                new_obj = new_obj2
            end

            if new_obj > best_obj
                best_obj = new_obj
                best_sol = deepcopy(temp_sol)
            end
        end
    end

    # for n in N
    #     if best_sol.planes[n].route[1] != best_sol.planes[n].route[end]
    #         println("HEY!!!")
    #     end
    # end

    return best_obj, best_sol
end

function ConstructLoop(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    N = data.N
    V = data.V
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        if m == 0
            base = Int(planes.planes[n].route[1])

            for v in V
                if v == base
                    continue
                end

                temp_sol = deepcopy(planes)
                Constructor(temp_sol.planes[n], v, 2, data)

                new_obj = obj(temp_sol, data, rt)

                temp_sol2 = deepcopy(temp_sol)
                UpdateTurnAroundTimes(temp_sol2, 1, maxTurnaround, data)
                new_obj2 = obj(temp_sol2, data, rt)

                if new_obj2 > new_obj
                    temp_sol = temp_sol2
                    new_obj = new_obj2
                end

                if new_obj > best_obj
                    best_obj = new_obj
                    best_sol = deepcopy(temp_sol)
                end
            end

            continue
        end

        for idx in 2:(m + 1)
            for v in V
                temp_sol = deepcopy(planes)

                if temp_sol.planes[n].route[idx] != v && temp_sol.planes[n].route[idx - 1] != v
                    Constructor(temp_sol.planes[n], v, idx, data)
                    new_obj = obj(temp_sol, data, rt)

                    temp_sol2 = deepcopy(temp_sol)
                    UpdateTurnAroundTimes(temp_sol2, idx, maxTurnaround, data)
                    new_obj2 = obj(temp_sol2, data, rt)

                    if new_obj < new_obj2
                        temp_sol = temp_sol2
                        new_obj = new_obj2
                    end

                    if new_obj > best_obj
                        best_obj = new_obj
                        best_sol = deepcopy(temp_sol)
                    end
                end
            end
        end
    end

    # for n in N
    #     if best_sol.planes[n].route[1] != best_sol.planes[n].route[end]
    #         println("HEY!!!")
    #     end
    # end

    return best_obj, best_sol
end

function Swap(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)
    N = data.N
    V = data.V
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Need at least 2 legs to have an internal stop to swap
        if m <= 1
            continue
        end

        # Internal route positions only (avoid route[1] base)
        for idx in 2:m
            # Neighbors around the original stop at idx
            left_v = planes.planes[n].route[idx - 1]
            old_v = planes.planes[n].route[idx]
            right_v = planes.planes[n].route[idx + 1]

            for v in V
                # Skip no-op and adjacent duplicates
                if v == old_v || v == left_v || v == right_v
                    continue
                end

                cand = deepcopy(planes)

                # 1) Insert new stop at idx
                Constructor(cand.planes[n], v, idx, data)

                # 2) Remove original stop (shifted to idx+1 after insertion)
                Destructor(cand.planes[n], [Int64(idx + 1)], data)

                new_obj = obj(cand, data, rt)

                # Candidate with updated turnaround times
                cand2 = deepcopy(cand)
                from_idx = Int64(max(1, idx - 1))
                UpdateTurnAroundTimes(cand2, from_idx, maxTurnaround, data)
                new_obj2 = obj(cand2, data, rt)

                if new_obj2 > new_obj
                    cand = cand2
                    new_obj = new_obj2
                end

                if new_obj > best_obj
                    best_obj = new_obj
                    best_sol = deepcopy(cand)
                end
            end
        end
    end

    # for n in N
    #     if best_sol.planes[n].route[1] != best_sol.planes[n].route[end]
    #         println("HEY!!!")
    #     end
    # end

    return best_obj, best_sol
end

function two_opt_Loop(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)
    N = data.N
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Need at least 3 legs to have two internal indices to reverse.
        if m <= 2
            continue
        end

        for i in 2:(m - 1)
            for j in (i + 1):m
                temp_sol = deepcopy(planes)
                two_opt(temp_sol.planes[n], i, j)

                temp_obj = obj(temp_sol, data, rt)

                # Re-evaluate with refreshed turnaround times from the changed area.
                temp_sol2 = deepcopy(temp_sol)
                from_idx = Int64(max(1, i - 1))
                UpdateTurnAroundTimes(temp_sol2, from_idx, maxTurnaround, data)
                temp_obj2 = obj(temp_sol2, data, rt)

                if temp_obj2 > temp_obj
                    temp_sol = temp_sol2
                    temp_obj = temp_obj2
                end

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                end
            end
        end
    end

    # for n in N
    #     if best_sol.planes[n].route[1] != best_sol.planes[n].route[end]
    #         println("HEY!!!")
    #     end
    # end
    
    return best_obj, best_sol
end



function RemoveVP(plane, data)

    ET = data.ET
    rt2 = data.rt
    te = data.te
    ec = data.ec
    e = data.e
    bmid = data.bmid
    bmin = data.bmin
    bmax = data.bmax
    

    
    m = Int(plane.flightLegs)

    for idx in 2:m
        temp_plane = deepcopy(plane)

        Destructor(temp_plane, [idx], data)

        if has_consecutive_duplicates(temp_plane.route)
            continue
        end

        Possible, SingleRouteBattery = Filter(temp_plane, rt2, e, te, ec, ET, bmin, bmid, bmax)

        if Possible 
            return temp_plane
        end
    end

    return nothing
end


function addVP(plane, data)


    ET = data.ET
    rt2 = data.rt
    te = data.te
    ec = data.ec
    e = data.e
    bmid = data.bmid
    bmin = data.bmin
    bmax = data.bmax
    V = data.V

    m = Int(plane.flightLegs)

    for idx in 2:(m+1)
        for v in V
            temp_plane = deepcopy(plane)
            if temp_plane.route[idx] != v && temp_plane.route[idx - 1] != v
                Constructor(temp_plane, v, idx, data)
                

                Possible, SingleRouteBattery = Filter(temp_plane, rt2, e, te, ec, ET, bmin, bmid, bmax)

                if Possible 
                    return temp_plane
                end

            end
        end
    end
    return nothing
end

function Repair(planes::allPlaneSolution, data, rt, newobj)
    ET = data.ET
    rt2 = data.rt
    te = data.te
    ec = data.ec
    e = data.e
    bmid = data.bmid
    bmin = data.bmin
    bmax = data.bmax
    N = data.N

    #Find planes with infeasible routes
   
    for (i, plane) in enumerate(planes.planes)
        Changed = true
        if plane.flightLegs > 1
            while newobj <= -900000 && Changed
                Possible, SingleRouteBattery = Filter(plane, rt2, e, te, ec, ET, bmin, bmid, bmax)

                if Possible
                    assign_Turnaround(planes, data) #Update to only do on infeasible routes.
                    # newobj = obj(planes, data, rt)
                    Changed = false
                else
                    plane_copy = deepcopy(plane)

                    if !SingleRouteBattery
                        newp = RemoveVP(plane, data)
                        if newp === nothing
                            newp = addVP(plane_copy, data)
                        end
                    else
                        newp = addVP(plane, data)
                        if newp === nothing
                            newp = RemoveVP(plane_copy, data)
                        end
                    end

                    if newp === nothing
                        return nothing
                    end

                    planes.planes[i] = newp
                    plane = newp
                end
            end
        end
    end

    return 0
end