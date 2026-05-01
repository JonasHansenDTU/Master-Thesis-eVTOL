
function weighted_choice(items::Vector{Int}, weights::Vector{Float64})
    #vertiports with higher demand weight are picked more often, while still keeping some randomness

    s = sum(weights)
    r = rand() * s
    cum = 0.0
    for (item, w) in zip(items, weights)
        cum += w
        if r <= cum
            return item
        end
    end
    return items[end]
end

function build_vertiport_weights(V, op, dp, A)
    #builds a vector of demand-based probabilities so that vertiports that appear often in passenger requests are more likely to be chosen 
    #when generating the initial chromosome

    counts = Dict(v => 1.0 for v in V)  # small base weight so every node is possible

    for a in A
        counts[op[a]] += 1.0
        counts[dp[a]] += 1.0
    end

    return [counts[v] for v in V]
end

function initial_chromosome_solution(data; maxLegs::Int=5, maxTurnaround::Int=30, debug_print::Bool=false)
    #creates a random initial route plan for each plane, starting and ending at its base,
    #with demand-biased intermediate vertiports and random flight legs and turnaround times, then packs everything into your chromosome structure.

    V = data.V
    N = data.N
    A = data.A
    bv = data.bv
    op = data.op
    dp = data.dp
    te = data.te
    dt = data.dt
    rt = data.rt
    w = data.w
    ET = data.ET

    weights = build_vertiport_weights(V, op, dp, A)

    planes = planeSolution[]
    nPlanes = length(N)

    choices = [i for i in 0:(maxLegs * nPlanes) if i != 1]
    total_flight_legs = rand(choices)

    allowed_legs = zeros(Int, nPlanes)

    while true
        remaining = total_flight_legs
        allowed_legs .= 0
        feasible_distribution = true

        for n in 1:(nPlanes - 1)
            feasible_choices = Int[]

            for x in 0:maxLegs
                rem_after = remaining - x
                if x != 1 && rem_after >= 0
                    max_possible_rest = maxLegs * (nPlanes - n)
                    if (rem_after == 0 || rem_after >= 2) && rem_after <= max_possible_rest
                        push!(feasible_choices, x)
                    end
                end
            end

            if isempty(feasible_choices)
                feasible_distribution = false
                break
            end

            choice = rand(feasible_choices)
            allowed_legs[n] = choice
            remaining -= choice
        end

        if feasible_distribution &&
           remaining <= maxLegs &&
           remaining != 1 &&
           remaining >= 0
            allowed_legs[nPlanes] = remaining
            allowed_legs = shuffle!(allowed_legs)
            break
        end
    end


    for n in 1:nPlanes
        base = bv[n]

        # choose number of legs
        flightLegs = allowed_legs[n]

        # route always starts at base
        route = Int32[base]

        if flightLegs > 0
            # choose intermediate nodes for first (flightLegs-1) legs
            current = base
            for k in 1:(flightLegs - 1)
                candidates = [v for v in V if v != current]
                cand_weights = [weights[findfirst(==(v), V)] for v in candidates]
                nxt = weighted_choice(candidates, cand_weights)
                push!(route, Int32(nxt))
                current = nxt
            end

            # final node forced back to base
            if current == base
                # if already at base, pick another node first when possible
                candidates = [v for v in V if v != base]
                if !isempty(candidates)
                    nxt = rand(candidates)
                    route[end] = Int32(nxt)
                end
            end
            push!(route, Int32(base))
        end

        turnaroundTime = Int32[]
        current_time = 0
        current_VP = base
        for k in 1:flightLegs
            common_a = [a for (a, v) in op if v == current_VP && dp[a] == route[k+1]]

            Candidate_Pass = [dt[a] for a in common_a if
                current_time + te - w <= dt[a] <= current_time + maxTurnaround &&
                dt[a] + rt[(op[a], dp[a])] <= ET
            ]
            if !isempty(Candidate_Pass)
                chosen_time = max(te,rand(Candidate_Pass) - current_time)
                push!(turnaroundTime, round(Int32, chosen_time))
            else
                # Choose rate (tune as needed)
                λ = 1.0 / (maxTurnaround - te)

                # Bounds as floats
                a = float(te)
                b = float(maxTurnaround)

                # Sample from truncated exponential using inverse CDF
                u = rand()
                x = a - log(1 - u * (1 - exp(-λ * (b - a)))) / λ

                # Convert to integer safely
                x_int = floor(Int, x)

                # Optional: clamp just in case of floating-point edge rounding
                x_int = clamp(x_int, Int(te), Int(maxTurnaround))

                # Store result
                push!(turnaroundTime, x_int)
            end
            current_time += turnaroundTime[end] + rt[(current_VP,route[(k+1)])]
            current_VP = route[(k+1)]
        end

        push!(planes, planeSolution(
            Int32(flightLegs),
            route,
            turnaroundTime
        ))
    end

    
    # ----- Check for optimal routing ------ #

    # Opt_evtol1 = [1,3,1]
    # Opt_evtol2 = [2, 3, 2]
    # # Opt_evtol3 = [5,2,3,5]

    # if planes[1].route == Opt_evtol1 && 
    #     planes[3].flightLegs == 0 &&
    #     planes[2].route == Opt_evtol2

    if debug_print
        # ----- Check for optimal routing ------ #

        Opt_evtol1 = [1,3,1]
        Opt_evtol2 = [2, 3, 2]
        # Opt_evtol3 = [5,2,3,5]

        if planes[1].route == Opt_evtol1 &&
            planes[3].flightLegs == 0 &&
            planes[2].route == Opt_evtol2

            rt = zeros(Int, Vmax, Vmax)
            for i in data.V, j in data.V
                rt[i, j] = data.rt[(i, j)]
            end

            UpdateTurnAroundTimes(allPlaneSolution(planes), 1, maxTurnaround, data)
            assignments, scheduled = assign_passengersV2(allPlaneSolution(planes), data, rt)
            Obj = obj(allPlaneSolution(planes), data, rt)
            println("Optimal Trips found, obj: $(Obj)")
            print_chromosome_table(allPlaneSolution(planes))
            print_assignments(assignments, data)
            print_schedule_pretty(scheduled)
        end

        #-----------------------------------------
        planes1 = allPlaneSolution(planes)
        for n in N
            if planes1.planes[n].route[1] != planes1.planes[n].route[end]
                println("HEY!!!")
            end
        end
    end
    #     rt = zeros(Int, Vmax, Vmax)
    #     for i in data.V, j in data.V
    #         rt[i, j] = data.rt[(i, j)]
    #     end

    #     UpdateTurnAroundTimes(allPlaneSolution(planes), 1, maxTurnaround, data)
    #     assignments, scheduled = assign_passengersV2(allPlaneSolution(planes), data, rt)
    #     Obj = obj(allPlaneSolution(planes), data, rt)
    #     println("Optimal Trips found, obj: $(Obj)")
    #     print_chromosome_table(allPlaneSolution(planes))
    #     print_assignments(assignments, data)
    #     print_schedule_pretty(scheduled) 
    # end


    return allPlaneSolution(planes)
end

function initial_chromosome_solution2(data, rt; maxLegs::Int=5, maxTurnaround::Int=30, top_c::Int=3, debug_print::Bool=false)
    V  = data.V
    N  = data.N
    A  = data.A
    bv = data.bv
    op = data.op
    dp = data.dp
    dt = data.dt
    q  = data.q
    so = data.so
    w  = Int(round(data.w))
    te = Int(round(data.te))
    fd = data.fd
    fs = data.fs
    cap_u = Int(round(data.cap_u))

    weights = build_vertiport_weights(V, op, dp, A)

    planes = Vector{Union{Nothing,planeSolution}}(undef, length(N))
    assignments = PassengerAssignment[]
    assigned_groups = Set{Int}()

    # ----------------------------
    # Random construction mode
    # ----------------------------
    early_first = rand() < 0.5

    passenger_value = Dict(a => (so[a] == 1 ? fs[(op[a], dp[a])] : fd[(op[a], dp[a])]) for a in A)

    if early_first
        sorted_groups = sort(A, by = a -> (dt[a], -passenger_value[a]))
    else
        sorted_groups = sort(A, by = a -> (-dt[a], -passenger_value[a]))
    end

    # ----------------------------
    # Random plane order
    # ----------------------------
    plane_order = shuffle(collect(1:length(N)))
    nPlanes = length(N)

    # ----------------------------
    # Distribute total flight legs
    # ----------------------------
    choices = [i for i in 0:(maxLegs * nPlanes) if i != 1]
    total_flight_legs = rand(choices)

    allowed_legs = zeros(Int, nPlanes)

    while true
        remaining = total_flight_legs
        allowed_legs .= 0
        feasible_distribution = true

        for n in 1:(nPlanes - 1)
            feasible_choices = Int[]

            for x in 0:maxLegs
                rem_after = remaining - x
                if x != 1 && rem_after >= 0
                    max_possible_rest = maxLegs * (nPlanes - n)
                    if (rem_after == 0 || rem_after >= 2) && rem_after <= max_possible_rest
                        push!(feasible_choices, x)
                    end
                end
            end

            if isempty(feasible_choices)
                feasible_distribution = false
                break
            end

            choice = rand(feasible_choices)
            allowed_legs[n] = choice
            remaining -= choice
        end

        if feasible_distribution &&
           remaining <= maxLegs &&
           remaining != 1 &&
           remaining >= 0
            allowed_legs[nPlanes] = remaining
            allowed_legs = shuffle!(allowed_legs)
            break
        end
    end

    # ----------------------------
    # Phase 1: build routes + seed assignments
    # ----------------------------
    for n in plane_order
        base = bv[n]
        flightLegs_target = allowed_legs[n]

        route = Int32[base]
        turnaroundTime = Int32[]

        current_vp = base
        current_time = 0
        legs_used = 0

        while legs_used < flightLegs_target
            candidates = NamedTuple[]

            for a in sorted_groups
                if a in assigned_groups
                    continue
                end

                # ------------------------
                # Direct option
                # ------------------------
                if current_vp == op[a]
                    if legs_used + 1 <= flightLegs_target
                        wait1 = max(te, Int(round(dt[a] - current_time)))
                        wait1 = min(wait1, maxTurnaround)
                        dep1 = current_time + wait1

                        if dt[a] <= dep1 <= dt[a] + w
                            score = passenger_value[a] - 0.1 * wait1
                            push!(candidates, (
                                score = score,
                                group = a,
                                mode = :direct,
                                midpoint = nothing,
                                reposition = false,
                                rep_wait = 0,
                                wait1 = wait1,
                                wait2 = 0
                            ))
                        end
                    end
                else
                    if legs_used + 2 <= flightLegs_target
                        rep_wait = te
                        arrival_origin = current_time + rep_wait + rt[current_vp, op[a]]
                        wait1 = max(te, Int(round(dt[a] - arrival_origin)))
                        wait1 = min(wait1, maxTurnaround)
                        dep1 = arrival_origin + wait1

                        if dt[a] <= dep1 <= dt[a] + w
                            score = passenger_value[a] - 0.1 * wait1 - 0.2 * rt[current_vp, op[a]]
                            push!(candidates, (
                                score = score,
                                group = a,
                                mode = :direct,
                                midpoint = nothing,
                                reposition = true,
                                rep_wait = rep_wait,
                                wait1 = wait1,
                                wait2 = 0
                            ))
                        end
                    end
                end

                # ------------------------
                # Stopover option
                # ------------------------
                if so[a] == 1
                    mids = [v for v in V if v != op[a] && v != dp[a]]
                    if !isempty(mids)
                        chosen_mid = weighted_choice(mids, [weights[findfirst(==(v), V)] for v in mids])

                        if current_vp == op[a]
                            if legs_used + 2 <= flightLegs_target
                                wait1 = max(te, Int(round(dt[a] - current_time)))
                                wait1 = min(wait1, maxTurnaround)
                                dep1 = current_time + wait1
                                arr1 = dep1 + rt[op[a], chosen_mid]
                                wait2 = te
                                dep2 = arr1 + wait2

                                if dt[a] <= dep1 <= dt[a] + w && dep2 >= arr1
                                    score = passenger_value[a] - 0.1 * (wait1 + wait2)
                                    push!(candidates, (
                                        score = score,
                                        group = a,
                                        mode = :stop,
                                        midpoint = chosen_mid,
                                        reposition = false,
                                        rep_wait = 0,
                                        wait1 = wait1,
                                        wait2 = wait2
                                    ))
                                end
                            end
                        else
                            if legs_used + 3 <= flightLegs_target
                                rep_wait = te
                                arrival_origin = current_time + rep_wait + rt[current_vp, op[a]]
                                wait1 = max(te, Int(round(dt[a] - arrival_origin)))
                                wait1 = min(wait1, maxTurnaround)
                                dep1 = arrival_origin + wait1
                                arr1 = dep1 + rt[op[a], chosen_mid]
                                wait2 = te
                                dep2 = arr1 + wait2

                                if dt[a] <= dep1 <= dt[a] + w && dep2 >= arr1
                                    score = passenger_value[a] - 0.1 * (wait1 + wait2) - 0.2 * rt[current_vp, op[a]]
                                    push!(candidates, (
                                        score = score,
                                        group = a,
                                        mode = :stop,
                                        midpoint = chosen_mid,
                                        reposition = true,
                                        rep_wait = rep_wait,
                                        wait1 = wait1,
                                        wait2 = wait2
                                    ))
                                end
                            end
                        end
                    end
                end
            end

            if isempty(candidates)
                break
            end

            candidates = sort(candidates, by = x -> x.score, rev = true)
            pool = candidates[1:min(top_c, length(candidates))]
            chosen = rand(pool)

            a = chosen.group

            if chosen.reposition
                push!(route, Int32(op[a]))
                push!(turnaroundTime, Int32(chosen.rep_wait))
                current_time += turnaroundTime[end] + rt[current_vp, op[a]]
                current_vp = op[a]
                legs_used += 1
            end

            if chosen.mode == :direct
                service_leg = legs_used + 1

                push!(route, Int32(dp[a]))
                push!(turnaroundTime, Int32(chosen.wait1))
                current_time += turnaroundTime[end] + rt[current_vp, dp[a]]
                current_vp = dp[a]
                legs_used += 1

                push!(assignments, PassengerAssignment(a, n, [service_leg]))
                push!(assigned_groups, a)

            elseif chosen.mode == :stop
                leg1 = legs_used + 1
                leg2 = legs_used + 2

                push!(route, Int32(chosen.midpoint))
                push!(turnaroundTime, Int32(chosen.wait1))
                current_time += turnaroundTime[end] + rt[current_vp, chosen.midpoint]
                current_vp = chosen.midpoint
                legs_used += 1

                push!(route, Int32(dp[a]))
                push!(turnaroundTime, Int32(chosen.wait2))
                current_time += turnaroundTime[end] + rt[current_vp, dp[a]]
                current_vp = dp[a]
                legs_used += 1

                push!(assignments, PassengerAssignment(a, n, [leg1, leg2]))
                push!(assigned_groups, a)
            end
        end

        # If no more passenger candidate: either fill spare legs randomly or return to base
        while legs_used < flightLegs_target - 1
            candidates = [v for v in V if v != current_vp]
            if isempty(candidates)
                break
            end

            nxt = weighted_choice(candidates, [weights[findfirst(==(v), V)] for v in candidates])
            push!(route, Int32(nxt))
            push!(turnaroundTime, Int32(rand(te:maxTurnaround)))
            current_time += turnaroundTime[end] + rt[current_vp, nxt]
            current_vp = nxt
            legs_used += 1
        end

        if current_vp != base
            push!(route, Int32(base))
            push!(turnaroundTime, Int32(rand(te:maxTurnaround)))
            current_time += turnaroundTime[end] + rt[current_vp, base]
            current_vp = base
            legs_used += 1
        end

        planes[n] = planeSolution(
            Int32(length(route) - 1),
            route,
            turnaroundTime
        )
    end

    evtols_init = allPlaneSolution(planeSolution[p for p in planes])
    scheduled = build_scheduled_legs(evtols_init, Int.(rt), cap_u)

    # ----------------------------
    # Phase 2: assign remaining groups on built schedule
    # ----------------------------
    assigned_groups = Set(ass.group for ass in assignments)

    remaining_capacity = Dict((leg.plane, leg.leg_index) => cap_u for leg in scheduled)
    exclusive_leg = Dict((leg.plane, leg.leg_index) => false for leg in scheduled)

    for ass in assignments
        a = ass.group
        for legidx in ass.legs
            remaining_capacity[(ass.plane, legidx)] -= q[a]
        end
        if so[a] == 0
            @assert length(ass.legs) == 1
            exclusive_leg[(ass.plane, ass.legs[1])] = true
        end
    end

    for a in sorted_groups
        if a in assigned_groups
            continue
        end

        earliest = Int(round(dt[a]))
        latest = earliest + w

        if so[a] == 0
            chosen = nothing
            best_arr = typemax(Int)

            for leg in scheduled
                key = (leg.plane, leg.leg_index)
                if leg.from == op[a] &&
                   leg.to == dp[a] &&
                   earliest <= leg.dep <= latest &&
                   remaining_capacity[key] == cap_u &&
                   !exclusive_leg[key]

                    if leg.arr < best_arr
                        best_arr = leg.arr
                        chosen = leg
                    end
                end
            end

            if chosen !== nothing
                key = (chosen.plane, chosen.leg_index)
                remaining_capacity[key] -= q[a]
                exclusive_leg[key] = true
                push!(assignments, PassengerAssignment(a, chosen.plane, [chosen.leg_index]))
                push!(assigned_groups, a)
            end
        else
            chosen_direct = nothing
            best_direct_arr = typemax(Int)

            for leg in scheduled
                key = (leg.plane, leg.leg_index)
                if leg.from == op[a] &&
                   leg.to == dp[a] &&
                   earliest <= leg.dep <= latest &&
                   remaining_capacity[key] >= q[a] &&
                   !exclusive_leg[key]

                    if leg.arr < best_direct_arr
                        best_direct_arr = leg.arr
                        chosen_direct = leg
                    end
                end
            end

            if chosen_direct !== nothing
                key = (chosen_direct.plane, chosen_direct.leg_index)
                remaining_capacity[key] -= q[a]
                push!(assignments, PassengerAssignment(a, chosen_direct.plane, [chosen_direct.leg_index]))
                push!(assigned_groups, a)
                continue
            end

            chosen_pair = nothing
            best_pair_arr = typemax(Int)

            for leg1 in scheduled
                key1 = (leg1.plane, leg1.leg_index)
                if leg1.from != op[a]
                    continue
                end
                if !(earliest <= leg1.dep <= latest)
                    continue
                end
                if remaining_capacity[key1] < q[a] || exclusive_leg[key1]
                    continue
                end

                for leg2 in scheduled
                    key2 = (leg2.plane, leg2.leg_index)
                    if leg2.plane != leg1.plane
                        continue
                    end
                    if leg2.leg_index != leg1.leg_index + 1
                        continue
                    end
                    if leg2.from != leg1.to
                        continue
                    end
                    if leg2.to != dp[a]
                        continue
                    end
                    if leg2.dep < leg1.arr
                        continue
                    end
                    if remaining_capacity[key2] < q[a] || exclusive_leg[key2]
                        continue
                    end

                    if leg2.arr < best_pair_arr
                        best_pair_arr = leg2.arr
                        chosen_pair = (leg1, leg2)
                    end
                end
            end

            if chosen_pair !== nothing
                leg1, leg2 = chosen_pair
                remaining_capacity[(leg1.plane, leg1.leg_index)] -= q[a]
                remaining_capacity[(leg2.plane, leg2.leg_index)] -= q[a]
                push!(assignments, PassengerAssignment(a, leg1.plane, [leg1.leg_index, leg2.leg_index]))
                push!(assigned_groups, a)
            end
        end
    end

    if debug_print
        for n in N
            if evtols_init.planes[n].route[1] != evtols_init.planes[n].route[end]
                println("HEY!!!")
            end
        end
    end

    return evtols_init, assignments, scheduled
end

function print_chromosome_table(evtols::allPlaneSolution)
    println("Chromosome table:")
    println("-----------------")

    for (n, plane) in enumerate(evtols.planes)
        print("eVTOL", n, ": ")
        print(plane.flightLegs, " | ")

        for v in plane.route
            print(v, " ")
        end

        print("| ")

        for t in plane.turnaroundTime
            print(t, " ")
        end

        println()
    end
end

function generate_best_initial_solutions(data, rt; n_runs::Int=1000, top_k::Int=10, top_c::Int=3,  maxLegs::Int=5, maxTurnaround::Int=30, print_each::Bool=false)
    results = NamedTuple[]

    for run in 1:n_runs
        evtols_init= initial_chromosome_solution(data; maxLegs=maxLegs, maxTurnaround=maxTurnaround)

        assignments, scheduled = assign_passengersV2(evtols_init, data, Int.(rt))

        P = FeasibilityCheck(Float32(data.bmax),Float32(data.bmid), Float32(data.bmin),data.dist,Float32(data.ec),Float32(data.battery_per_km),
                    evtols_init,Int.(rt),Int(round(data.ET)),maximum(Int.(data.T)),maximum(data.V),data.cap_v, data.b_penalty)

        fitness = fitnessFunction(evtols_init,assignments,Float32(data.bmax), Float32(data.bmid), Float32(data.bmin),data.dist, Float32(data.ec),
                        Float32(data.battery_per_km), Int.(rt), Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V), data.cap_v, data)

        push!(results, (run = run, fitness = fitness, evtols = evtols_init, assignments = assignments, scheduled = scheduled, P = P))

        if print_each
            println("Run $run | fitness = $fitness | P = $P")
        end
    end

    # Sort by fitness descending, since higher fitness is better
    sorted_results = sort(results, by = x -> x.fitness, rev = true)

    # Return only the best top_k results
    return sorted_results[1:min(top_k, length(sorted_results))]
end
