mutable struct SingleRoutePoolEntry
    plane::planeSolution
    score::Float64
    source_plane::Int64
end

function single_route_key(plane::planeSolution)
    return (Int(plane.flightLegs), Tuple(Int.(plane.route)))
end

function route_contains_sequence(route, seq::Vector{Int})
    r = Int.(route)
    n = length(r)
    m = length(seq)
    if n < m
        return false
    end
    for i in 1:(n - m + 1)
        match = true
        for j in 1:m
            if r[i + j - 1] != seq[j]
                match = false
                break
            end
        end
        if match
            return true
        end
    end
    return false
end

function score_single_plane_solution(plane::planeSolution, data, rt)
    single_sol = allPlaneSolution([deepcopy(plane)])
    return obj(single_sol, data, rt)
end

function is_single_plane_feasible(plane::planeSolution, data, rt)
    single_sol = allPlaneSolution([deepcopy(plane)])
    P = FeasibilityCheck(
        Float32(data.bmax),
        Float32(data.bmid),
        Float32(data.bmin),
        data.dist,
        Float32(data.ec),
        Float32(data.battery_per_km),
        single_sol,
        Int.(rt),
        Int(round(data.ET)),
        Int(maximum(data.T)),
        Int(maximum(data.V)),
        data.cap_v,
        data.b_penalty
    )
    return sum(P) < 1000000
end

function check_turnaround_times(sol::allPlaneSolution, sol_name::String="")
    """Check if any plane in sol has turnaround times < 30 min at index >= 2."""
    issues = false
    label = isempty(sol_name) ? "" : "[$sol_name] "
    
    for (pidx, plane) in enumerate(sol.planes)
        if plane.flightLegs <= 0
            continue
        end
        
        for tur_idx in 2:length(plane.turnaroundTime)
            if plane.turnaroundTime[tur_idx] < 30
                println("WARNING: $(label)Plane $pidx has turnaround time < 30 min at index $tur_idx: $(plane.turnaroundTime[tur_idx]) min")
                issues = true
            end
        end
    end
    
    return issues
end

function collect_feasible_single_plane_routes(sol::allPlaneSolution, data, rt; pool=SingleRoutePoolEntry[], max_size::Int=50, max_duplicates::Int=3)
    # Store up to max_duplicates entries per route shape (key)
    routes_by_key = Dict{Any, Vector{SingleRoutePoolEntry}}()

    # Add existing pool entries
    for entry in pool
        key = single_route_key(entry.plane)
        if !haskey(routes_by_key, key)
            routes_by_key[key] = SingleRoutePoolEntry[]
        end
        push!(routes_by_key[key], deepcopy(entry))
    end

    # Add new entries from current solution
    for (idx, plane) in enumerate(sol.planes)
        if plane.flightLegs <= 0
            continue
        end

        if !is_single_plane_feasible(plane, data, rt)
            continue
        end

        key = single_route_key(plane)
        score = score_single_plane_solution(plane, data, rt)
        new_entry = SingleRoutePoolEntry(deepcopy(plane), score, Int64(idx))

        if !haskey(routes_by_key, key)
            routes_by_key[key] = SingleRoutePoolEntry[]
        end
        push!(routes_by_key[key], new_entry)
    end

    # For each route shape, keep only top max_duplicates by score
    for key in keys(routes_by_key)
        entries = routes_by_key[key]
        sort!(entries, by = e -> e.score, rev = true)
        if length(entries) > max_duplicates
            routes_by_key[key] = entries[1:max_duplicates]
        end
    end

    # Flatten to single vector and sort by score
    pool = SingleRoutePoolEntry[]
    for entries in values(routes_by_key)
        append!(pool, entries)
    end

    sort!(pool, by = entry -> entry.score, rev = true)

    if length(pool) > max_size
        resize!(pool, max_size)
    end

    return pool
end

function build_pool_candidate(pool::Vector{SingleRoutePoolEntry}, data, rt; max_routes::Int=4)
    if isempty(pool) || isempty(data.N)
        return nothing
    end

    # empty_plane placeholder not used; create per-slot with base port when needed

    # Build a fresh candidate from the pool only. Use data.N for plane indices
    n_planes = length(data.N)
    cand_planes = Vector{planeSolution}(undef, n_planes)

    # Track which route-keys we've already applied so we don't pick the same route twice.
    selected_keys = Set{Any}()

    # Process planes in random order for more diversity (use positions)
    plane_indices = collect(1:n_planes)
    shuffle!(plane_indices)

    for pos in plane_indices
        # map position -> plane id
        pid = data.N[pos]
        # Get baseport of current plane from data.bv (base vertiport mapping)
        base_port = Int(data.bv[pid])

        # Find all routes in pool that start at the same baseport
        compatible_routes = [entry for entry in pool if Int(entry.plane.route[1]) == base_port]

        if isempty(compatible_routes)
            cand_planes[pos] = planeSolution(Int32(0), Int32[base_port], Int32[])
            continue
        end

        # Preserve the already-sorted pool order; filtering keeps that order.
        filtered = [e for e in compatible_routes if !(single_route_key(e.plane) in selected_keys)]
        if isempty(filtered)
            cand_planes[pos] = planeSolution(Int32(0), Int32[base_port], Int32[])
            continue
        end

        # Randomly pick from the top 10 (or fewer if not available)
        # top_k = min(100, length(filtered))
        choice = filtered

        # Apply chosen route and mark its key as used
        key = single_route_key(choice.plane)
        push!(selected_keys, key)
        cand_planes[pos] = deepcopy(choice.plane)
    end

    cand = allPlaneSolution(cand_planes)

    # Check if this candidate matches the target solution pattern
    # target_keys = [
    #     (2, (1, 3, 1)),
    #     (2, (2, 3, 2)),
    #     (2, (5, 3, 5))
    # ]
    
    # cand_keys = [single_route_key(plane) for plane in cand.planes]
    # if cand_keys == target_keys
    #     println("\n========== TARGET SOLUTION GENERATED! ==========")
    #     println("Found target solution from pool candidate!")
    #     for (i, plane) in enumerate(cand.planes)
    #         if plane.flightLegs > 0
    #             println("  eVTOL$(i): $(plane.flightLegs) | $(Int.(plane.route)) | $(Int.(plane.turnaroundTime))")
    #         else
    #             println("  eVTOL$(i): $(plane.flightLegs) | $(Int.(plane.route)) |")
    #         end
    #     end
    #     try
    #         target_obj = obj(cand, data, Int.(rt))
    #         println("\nObjective: $target_obj")
    #     catch e
    #         println("Error computing objective: $e")
    #     end
        
    #     println("\nFeasibility Check:")
    #     try
    #         P = FeasibilityCheck(
    #             Float32(data.bmax),
    #             Float32(data.bmid),
    #             Float32(data.bmin),
    #             data.dist,
    #             Float32(data.ec),
    #             Float32(data.battery_per_km),
    #             cand,
    #             Int.(rt),
    #             Int(round(data.ET)),
    #             Int(maximum(data.T)),
    #             Int(maximum(data.V)),
    #             data.cap_v,
    #             data.b_penalty
    #         )
    #         println("  Penalty vector P: $P")
    #         println("  Battery: $(P[1]), CompletionTime: $(P[2]), VertiportCapacity: $(P[3])")
    #     catch e
    #         println("Error computing feasibility check: $e")
    #     end
        
    #     println("\nPassenger Assignment:")
    #     try
    #         assignments, scheduled = assign_passengersV2(cand, data, Int.(rt))
    #         print_assignments(assignments, data)
    #     catch e
    #         println("Error computing passenger assignment: $e")
    #     end
    #     println("=============================================\n")
    # end

    return cand
end

function Heuristic(maxTurnaround::Int64, MaxTime::Int32, data, rt, top_c)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    candiateroutes = Candidate_Route(data)

    
    while elapsed <= Float64(MaxTime)
        nr = 1
        idx = rand(1:nr)
        T = 100
        Best_sols = generate_best_initial_solutions(data, rt, candiateroutes; n_runs = 20, top_k = nr, top_c, maxLegs=6, maxTurnaround)

        temp_obj = Best_sols[idx].fitness
        temp_sol = deepcopy(Best_sols[idx].evtols)
        if !check_turnaround_times(temp_sol, "Initial solution")
            println("HEY!")
        end

        if temp_obj > best_obj
            best_obj = temp_obj
            best_sol = deepcopy(temp_sol)
            println("New best Obj $(best_obj)")
            println("Method used: Initial Heuristic")
        end

        improvement = true
        while improvement && elapsed <= Float64(MaxTime)
            improvement = false
            method_used = 0

            des_obj, des_sol = DestructLoop(temp_sol, maxTurnaround, temp_obj, data, rt)
            con_obj, con_sol = ConstructLoop(temp_sol, maxTurnaround, temp_obj, data, rt)
            swap_obj, swap_sol = Swap(temp_sol, maxTurnaround, temp_obj, data, rt)
            two_opt_obj, two_opt_sol = two_opt_Loop(temp_sol, maxTurnaround, temp_obj, data, rt)


            ### ------ Repair functionTest ------------

            if des_obj <= -100000
                Repair(des_sol, data, rt, des_obj)
            end        
            if con_obj <= -100000
                Repair(con_sol, data, rt, con_obj)
            end     
            if swap_obj <= -100000
                Repair(swap_sol, data, rt, swap_obj)
            end     
            if two_opt_obj <= -100000
                Repair(two_opt_sol, data, rt, two_opt_obj)
            end     

            ### ---------------------------------------

            cand_obj = des_obj
            cand_sol = des_sol

            if con_obj > cand_obj || rand() < exp(-1/(T))
                cand_obj = con_obj
                cand_sol = con_sol
                method_used = 1
            end

            if swap_obj > cand_obj || rand() < exp(-1/(T))
                cand_obj = swap_obj
                cand_sol = swap_sol
                method_used = 2
            end

            if two_opt_obj > cand_obj || rand() < exp(-1/(T))
                cand_obj = two_opt_obj
                cand_sol = two_opt_sol
                method_used = 3
            end

            if cand_obj > temp_obj || rand() < exp(-1/(T))
                temp_obj = cand_obj
                temp_sol = deepcopy(cand_sol)
                if !check_turnaround_times(temp_sol, "After neighborhood")
                    println("HEY!")
                end
                improvement = true

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                    println("New best Obj $(best_obj)")
                    println("Method used: $(("Destructor", "Constructor", "Swap", "2-opt")[method_used + 1])")
                end
            end
            T = T * 0.95

            elapsed = (time_ns() - start_ns) / 1e9
        end

        iterations += 1
        if iterations % 100 == 0
            println("iteration: $(iterations)")
            println("Obj In this iteration $(temp_obj)")
        end
        elapsed = (time_ns() - start_ns) / 1e9
    end

    return best_obj, best_sol, iterations
end

function HeuristicSA(maxTurnaround::Int64, MaxTime::Int32, data, rt, top_c)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    candiateroutes = Candidate_Route(data)
    single_route_pool = SingleRoutePoolEntry[]


    max_size = length(data.N)*length(data.V) 

    
    while elapsed <= Float64(MaxTime)
        nr = 5
        idx = rand(1:nr)
        T = 50
        Best_sols = generate_best_initial_solutions(data, rt, candiateroutes; n_runs = 10, top_k = nr, top_c, maxLegs=6, maxTurnaround)

        temp_obj = Best_sols[idx].fitness
        temp_sol = deepcopy(Best_sols[idx].evtols)

        for best_sol_candidate in Best_sols
            single_route_pool = collect_feasible_single_plane_routes(best_sol_candidate.evtols, data, rt; pool=single_route_pool, max_size=max_size)
        end

        if temp_obj > best_obj
            best_obj = temp_obj
            best_sol = deepcopy(temp_sol)
            println("New best Obj $(best_obj)")
            println("Method used: Initial Heuristic")
        end

        single_route_pool = collect_feasible_single_plane_routes(temp_sol, data, rt; pool=single_route_pool, max_size=max_size)

        if !isempty(single_route_pool)
            pool_sol = build_pool_candidate(single_route_pool, data, rt; max_routes=min(2, length(single_route_pool)))
            if pool_sol !== nothing
                pool_obj = obj(pool_sol, data, rt)
                if pool_obj > temp_obj || rand() < 0.1
                    temp_sol = pool_sol
                    temp_obj = pool_obj
                    single_route_pool = collect_feasible_single_plane_routes(temp_sol, data, rt; pool=single_route_pool, max_size=max_size)
                end
            end
        end

        improvement = true
        while improvement && elapsed <= Float64(MaxTime)
            improvement = false
            method_used = 0

            cand_obj = temp_obj
            cand_sol = deepcopy(temp_sol)
            n_perm = 2

            for _ in 1:n_perm
                if rand(Bool)
                    cand_obj, cand_sol = DestructSA(cand_sol, maxTurnaround, cand_obj, data, rt)
                    next_method = 0
                else
                    cand_obj, cand_sol = ConstructSA(cand_sol, maxTurnaround, cand_obj, data, rt)
                    next_method = 1
                end
            end

            # out = 0
            
            # if cand_obj <= -900000
                    
            #     out = Repair(cand_sol, data, rt, cand_obj)
            #     if isnothing(out)
            #         continue
            #     end
            #     cand_obj = obj(cand_sol, data, rt)

            # end 


            if cand_obj > temp_obj || (rand() < exp((temp_obj - cand_obj) / T) && temp_obj > cand_obj)
                temp_obj = cand_obj
                temp_sol = deepcopy(cand_sol)
                improvement = true
                single_route_pool = collect_feasible_single_plane_routes(temp_sol, data, rt; pool=single_route_pool, max_size=max_size)

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                    println("New best Obj $(best_obj)")
                    println("Method used: $(("Destructor", "Constructor")[method_used + 1])")
                end
            end

            elapsed = (time_ns() - start_ns) / 1e9
            T = max(T *0.9, 0.0001)
        end

        iterations += 1
        if iterations % 100 == 0
            println("iteration: $(iterations)")
            println("Obj In this iteration $(temp_obj)")
        end
        elapsed = (time_ns() - start_ns) / 1e9
    end

    println("\n========== FINAL POOL SOLUTIONS ==========")
    println("Total pool entries: $(length(single_route_pool))")
    for (idx, entry) in enumerate(single_route_pool)
        key = single_route_key(entry.plane)
        println("  [$idx] flightLegs=$(key[1]), route=$(key[2]), score=$(entry.score)")
    end
    println("=========================================\n")

    return best_obj, best_sol, iterations
end




