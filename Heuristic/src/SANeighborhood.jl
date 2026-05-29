function DestructSA(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt, n_trials::Int=1)

    N = data.N
    trials = max(1, n_trials)

    best_obj = init_obj
    best_sol = planes
    found_candidate = false

    for _ in 1:trials
        n = rand(N)
        m = Int(planes.planes[n].flightLegs)

        if m == 0
            continue
        end

        if m == 1
            temp_sol = deepcopy(planes)
            base = temp_sol.planes[n].route[1]
            resize!(temp_sol.planes[n].route, 1)
            temp_sol.planes[n].route[1] = base
            empty!(temp_sol.planes[n].turnaroundTime)
            temp_sol.planes[n].flightLegs = Int32(0)

            new_obj = obj(temp_sol, data, rt)

            if !found_candidate || new_obj > best_obj
                best_obj = new_obj
                best_sol = temp_sol
                found_candidate = true
            end
            continue
        end

        # Do not delete first node (route[1]); last route node is not in this range anyway
        idx = rand(2:m)
        temp_sol = deepcopy(planes)

        Destructor(temp_sol.planes[n], [idx], data)

        if has_consecutive_duplicates(temp_sol.planes[n].route)
            continue
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

        if !found_candidate || new_obj > best_obj
            best_obj = new_obj
            best_sol = temp_sol
            found_candidate = true
        end
    end

    if found_candidate
        return best_obj, best_sol
    end

    return init_obj, planes
end

function ConstructSA(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    N = data.N
    V = data.V
    best_obj = init_obj
    best_sol = deepcopy(planes)

    n = rand(N)
    m = Int(planes.planes[n].flightLegs)

    if m == 0
        base = Int(planes.planes[n].route[1])
        end_VP_choices = data.end_vp[base]
        
        # Filter out base to prevent single-leg base-to-base
        end_candidates = [v for v in end_VP_choices if v != base]
        
        if isempty(end_candidates)
            # Fallback: build 2-leg route base -> random_non_base -> base
            non_base = [v for v in V if v != base]
            if isempty(non_base)
                return init_obj, planes
            end

            temp_sol = deepcopy(planes)
            te32 = Int32(data.te)
            mid = Int32(rand(non_base))

            insert!(temp_sol.planes[n].route, 2, mid)
            insert!(temp_sol.planes[n].turnaroundTime, 1, te32)
            insert!(temp_sol.planes[n].route, 3, Int32(base))
            insert!(temp_sol.planes[n].turnaroundTime, 2, te32)
            temp_sol.planes[n].flightLegs = Int32(2)

            new_obj = obj(temp_sol, data, rt)

            temp_sol2 = deepcopy(temp_sol)
            UpdateTurnAroundTimes(temp_sol2, 1, maxTurnaround, data)
            new_obj2 = obj(temp_sol2, data, rt)

            if new_obj2 > new_obj
                temp_sol = temp_sol2
                new_obj = new_obj2
            end

            return new_obj, temp_sol
        end
        endpoint = rand(end_candidates)

        temp_sol = deepcopy(planes)
        
        # Build single-leg route: [base, endpoint]
        te32 = Int32(data.te)
        insert!(temp_sol.planes[n].route, 2, Int32(endpoint))
        insert!(temp_sol.planes[n].turnaroundTime, 1, te32)
        temp_sol.planes[n].flightLegs = 1

        new_obj = obj(temp_sol, data, rt)

        temp_sol2 = deepcopy(temp_sol)
        UpdateTurnAroundTimes(temp_sol2, 1, maxTurnaround, data)
        new_obj2 = obj(temp_sol2, data, rt)

        if new_obj2 > new_obj
            temp_sol = temp_sol2
            new_obj = new_obj2
        end


        return new_obj, temp_sol
    end

    idx = rand(2:(m+1))
    
    # For final position insertions, pick from allowed endpoints
    if idx == m + 1
        base = Int(planes.planes[n].route[1])
        end_VP_choices = data.end_vp[base]
        current_last = Int(planes.planes[n].route[m])
        
        # Filter out current last node to prevent duplicates
        v_candidates = [v for v in end_VP_choices if v != current_last]
        
        if isempty(v_candidates)
            return init_obj, planes
        end
        
        v = rand(v_candidates)
    else
        v = rand(V)
    end
    
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

        return new_obj, temp_sol

    else
        return init_obj, planes 
    end


    # for n in N
    #     if best_sol.planes[n].route[1] != best_sol.planes[n].route[end]
    #         println("HEY!!!")
    #     end
    # end

end



function SegmentExchange(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    # Need at least two planes for a cross-route exchange.
    if length(planes.planes) <= 1
        return init_obj, planes
    end

    # A segment is chosen from internal route nodes only, so each route must have
    # at least one internal node (flightLegs >= 2).
    eligible = [i for i in eachindex(planes.planes) if Int(planes.planes[i].flightLegs) >= 2]
    if length(eligible) < 2
        return init_obj, planes
    end

    n1 = rand(eligible)
    n2_choices = [i for i in eligible if i != n1]
    n2 = rand(n2_choices)

    m1 = Int(planes.planes[n1].flightLegs)
    m2 = Int(planes.planes[n2].flightLegs)

    # Random contiguous segment in internal nodes: [2, m].
    s1 = rand(2:m1)
    e1 = rand(s1:m1)
    s2 = rand(2:m2)
    e2 = rand(s2:m2)

    temp_sol = deepcopy(planes)

    seg1 = copy(temp_sol.planes[n1].route[s1:e1])
    seg2 = copy(temp_sol.planes[n2].route[s2:e2])

    splice!(temp_sol.planes[n1].route, s1:e1, seg2)
    splice!(temp_sol.planes[n2].route, s2:e2, seg1)

    # Reject moves that create immediate repeated vertiports.
    if has_consecutive_duplicates(temp_sol.planes[n1].route) ||
       has_consecutive_duplicates(temp_sol.planes[n2].route)
        return init_obj, planes
    end

    temp_sol.planes[n1].flightLegs = Int32(length(temp_sol.planes[n1].route) - 1)
    temp_sol.planes[n2].flightLegs = Int32(length(temp_sol.planes[n2].route) - 1)

    # Rebuild turnaround vectors to match potentially changed route lengths.
    empty!(temp_sol.planes[n1].turnaroundTime)
    for _ in 1:Int(temp_sol.planes[n1].flightLegs)
        push!(temp_sol.planes[n1].turnaroundTime, Int32(data.te))
    end

    empty!(temp_sol.planes[n2].turnaroundTime)
    for _ in 1:Int(temp_sol.planes[n2].flightLegs)
        push!(temp_sol.planes[n2].turnaroundTime, Int32(data.te))
    end

    new_obj = obj(temp_sol, data, rt)

    temp_sol2 = deepcopy(temp_sol)
    from_idx = Int64(max(1, min(s1, s2) - 1))
    UpdateTurnAroundTimes(temp_sol2, from_idx, maxTurnaround, data)
    new_obj2 = obj(temp_sol2, data, rt)

    if new_obj2 > new_obj
        temp_sol = temp_sol2
        new_obj = new_obj2
    end

    return new_obj, temp_sol
end