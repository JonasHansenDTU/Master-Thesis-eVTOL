mutable struct SingleRoutePoolEntry
    plane::planeSolution
    score::Float64
    source_plane::Int64
    assignments::Vector{Any}
    diversity_score::Float64
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
    
    # Get passenger assignments for this single plane route
    assignments, scheduled = assign_passengersV2(single_sol, data, Int.(rt))
    
    # Compute the objective value directly
    fitness = fitnessFunction(
        single_sol,
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
    
    # Return both the fitness value and the passenger assignments
    return fitness, assignments
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

function check_route_endpoints(sol::allPlaneSolution, data, sol_name::String="")
    issues = false
    label = isempty(sol_name) ? "" : "[$sol_name] "

    for (pidx, plane) in enumerate(sol.planes)
        if plane.flightLegs <= 0
            continue
        end

        base = Int(plane.route[1])
        endpoint = Int(plane.route[end])
        end_VP_choices = get(data.end_vp, base, Int[])

        if !(endpoint in end_VP_choices)
            println("WARNING: $(label)Plane $pidx ends at vertiport $endpoint, not in end_vp[$base]=$(sort(collect(end_VP_choices)))")
            println("         Route: $(Int.(plane.route))")
            issues = true
        end
    end

    return issues
end

function get_passenger_set(assignments::Vector{Any})
    """Extract set of passenger group IDs from assignments."""
    passenger_groups = Set{Int}()
    for ass in assignments
        push!(passenger_groups, ass.group)
    end
    return passenger_groups
end

function get_served_od_pairs(assignments::Vector{Any}, data)
    """Extract set of (origin, destination) tuples served by assignments."""
    od_pairs = Set{Tuple{Int,Int}}()
    op = data.op
    dp = data.dp
    for ass in assignments
        group_id = ass.group
        origin = op[group_id]
        destination = dp[group_id]
        push!(od_pairs, (origin, destination))
    end
    return od_pairs
end

function compute_diversity_metrics(entry::SingleRoutePoolEntry, pool::Vector{SingleRoutePoolEntry}, data)
    """Compute hybrid diversity score: coverage_spread + od_region_diversity.
    
    Returns a score in [0, 1] where higher is more diverse.
    Diversity is computed relative to OTHER ROUTES FROM THE SAME BASE vertiport.
    - coverage_spread: How many unique passengers NOT served by other same-base routes
    - od_region_diversity: How many unique OD pairs NOT served by other same-base routes
    
    Metric: diversity_score = (unique_passengers / total_passengers) * 0.5 + (unique_od_pairs / total_od_pairs) * 0.5
    
    This prevents routes from the same base from all targeting the same passengers.
    """
    if isempty(entry.assignments)
        return 0.0  # Empty routes have no diversity benefit
    end
    
    # Get this entry's passengers and OD pairs
    entry_passengers = get_passenger_set(entry.assignments)
    entry_od_pairs = get_served_od_pairs(entry.assignments, data)
    
    # Get this entry's base vertiport
    entry_base = Int(entry.plane.route[1])
    
    # Get all passengers and OD pairs from OTHER routes with the SAME base vertiport
    # This ensures diversity is computed relative to competing routes from the same base
    pool_passengers = Set{Int}()
    pool_od_pairs = Set{Tuple{Int,Int}}()
    for other_entry in pool
        if other_entry === entry  # Skip self
            continue
        end
        # Only compare against routes from the same base vertiport
        if Int(other_entry.plane.flightLegs) > 0
            other_base = Int(other_entry.plane.route[1])
            if other_base == entry_base
                union!(pool_passengers, get_passenger_set(other_entry.assignments))
                union!(pool_od_pairs, get_served_od_pairs(other_entry.assignments, data))
            end
        end
    end
    
    # Compute coverage spread: unique passengers not in pool
    unique_passengers = setdiff(entry_passengers, pool_passengers)
    coverage_spread = isempty(entry_passengers) ? 0.0 : length(unique_passengers) / length(entry_passengers)
    
    # Compute OD region diversity: unique OD pairs not in pool
    unique_od_pairs = setdiff(entry_od_pairs, pool_od_pairs)
    od_diversity = isempty(entry_od_pairs) ? 0.0 : length(unique_od_pairs) / length(entry_od_pairs)
    
    # Hybrid score: weighted average
    diversity_score = 0.5 * coverage_spread + 0.5 * od_diversity
    
    return diversity_score
end

function collect_feasible_single_plane_routes(sol::allPlaneSolution, data, rt; pool=SingleRoutePoolEntry[], max_size::Int=50, max_duplicates::Int=3, debug::Bool=false)
    # Store up to max_duplicates entries per route shape (key)
    routes_by_key = Dict{Any, Vector{SingleRoutePoolEntry}}()

    maxTurnaround = Int64(data.ET)
    
    # if debug
    #     active_routes = sum((1 for p in sol.planes if Int(p.flightLegs) > 0); init=0)
    #     println("\n--- Processing solution with $active_routes active routes (pool has $(length(pool)) entries) ---")
    # end

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
        fitness, assignments = score_single_plane_solution(plane, data, rt)


        plane_sol = allPlaneSolution([plane])

        best_obj = fitness
        best_sol = deepcopy(plane_sol)

        improved = true
        while improved
            improved = false
            cand_obj, cand_sol  = DestructLoop(plane_sol, maxTurnaround, fitness, data, rt)
            con_obj, con_sol  = ConstructLoop(plane_sol, maxTurnaround, fitness, data, rt)

            if con_obj > cand_obj
                cand_obj = con_obj
                cand_sol = con_sol
            end

            if cand_obj > best_obj
                best_obj = cand_obj
                best_sol = deepcopy(cand_sol)
                improved = true
            end

        end

        if best_obj > fitness
            plane = deepcopy(best_sol.planes[1])
        end

        diversity = compute_diversity_metrics(SingleRoutePoolEntry(deepcopy(plane), fitness, Int64(idx), deepcopy(assignments), 0.0), pool, data)
        new_entry = SingleRoutePoolEntry(deepcopy(plane), fitness, Int64(idx), deepcopy(assignments), diversity)
        
        # Debug output: show what routes are being added to the pool
        if debug && Int(plane.flightLegs) > 0
            base_port = Int(plane.route[1])
            passenger_groups = [ass.group for ass in assignments]
            groups_str = isempty(passenger_groups) ? "none" : join(passenger_groups, ",")
            # println("  [ADD] Plane $idx (base $base_port): legs=$(Int(plane.flightLegs)), route=$(Int.(plane.route)), diversity=$(round(diversity*100; digits=1))%, passengers=[$groups_str]")
        end

        if !haskey(routes_by_key, key)
            routes_by_key[key] = SingleRoutePoolEntry[]
        end
        push!(routes_by_key[key], new_entry)
    end

    # Keep an explicit empty route for each base so inactive planes are always possible.
    unique_bases = unique(Int(data.bv[n]) for n in data.N)
    for base in unique_bases
        empty_plane = planeSolution(Int32(0), Int32[base], Int32[])
        empty_key = single_route_key(empty_plane)

        if !haskey(routes_by_key, empty_key)
            routes_by_key[empty_key] = SingleRoutePoolEntry[]
        end

        if isempty(routes_by_key[empty_key])
            push!(routes_by_key[empty_key], SingleRoutePoolEntry(empty_plane, 0.0, Int64(0), [], 0.0))
        else
            routes_by_key[empty_key][1] = SingleRoutePoolEntry(empty_plane, 0.0, Int64(0), [], 0.0)
            resize!(routes_by_key[empty_key], 1)
        end
    end

    # For each route shape, keep only top max_duplicates by diversity first, then fitness.
    # Empty routes are kept as exactly one entry per base.
    # Ranking: (diversity_score DESC, fitness_score DESC)
    for key in keys(routes_by_key)
        entries = routes_by_key[key]
        # Sort by (-diversity_score, -fitness_score) to sort both descending
        sort!(entries, by = e -> (-e.diversity_score, -e.score))

        if key[1] == 0
            routes_by_key[key] = entries[1:1]
        elseif length(entries) > max_duplicates
            routes_by_key[key] = entries[1:max_duplicates]
        end
    end

    # Flatten to single vector and split by active/inactive routes.
    all_entries = SingleRoutePoolEntry[]
    for entries in values(routes_by_key)
        append!(all_entries, entries)
    end

    empty_entries = [e for e in all_entries if Int(e.plane.flightLegs) == 0]
    sort!(empty_entries, by = e -> Int(e.plane.route[1]))

    non_empty_entries = [e for e in all_entries if Int(e.plane.flightLegs) > 0]
    # Maintain diversity-first ranking: sort by (-diversity, -fitness) to preserve priority
    sort!(non_empty_entries, by = e -> (-e.diversity_score, -e.score))

    # Group active routes by start base to preserve cross-base diversity.
    entries_by_base = Dict{Int, Vector{SingleRoutePoolEntry}}()
    for entry in non_empty_entries
        base_port = Int(entry.plane.route[1])
        if !haskey(entries_by_base, base_port)
            entries_by_base[base_port] = SingleRoutePoolEntry[]
        end
        push!(entries_by_base[base_port], entry)
    end

    for entries in values(entries_by_base)
        # Maintain diversity-first ranking within each base
        sort!(entries, by = e -> (-e.diversity_score, -e.score))
    end

    # Reserve room for empties, then allocate active routes with a soft per-base quota.
    keep = SingleRoutePoolEntry[]
    append!(keep, empty_entries)

    remaining_capacity = max_size - length(keep)
    if remaining_capacity <= 0
        if length(keep) > max_size
            resize!(keep, max_size)
        end
        return keep
    end

    non_empty_bases = sort(collect(keys(entries_by_base)))
    nbases = length(non_empty_bases)

    selected_non_empty = SingleRoutePoolEntry[]
    remainder_non_empty = SingleRoutePoolEntry[]

    if nbases > 0
        per_base_quota = max(1, Int(floor(remaining_capacity / nbases)))

        for base in non_empty_bases
            entries = entries_by_base[base]
            take_n = min(per_base_quota, length(entries))
            if take_n > 0
                append!(selected_non_empty, entries[1:take_n])
            end
            if length(entries) > take_n
                append!(remainder_non_empty, entries[take_n+1:end])
            end
        end

        sort!(selected_non_empty, by = e -> e.score, rev = true)
        if length(selected_non_empty) > remaining_capacity
            resize!(selected_non_empty, remaining_capacity)
            append!(keep, selected_non_empty)
            return keep
        end

        append!(keep, selected_non_empty)

        extra_capacity = max_size - length(keep)
        if extra_capacity > 0 && !isempty(remainder_non_empty)
            sort!(remainder_non_empty, by = e -> e.score, rev = true)
            append!(keep, remainder_non_empty[1:min(extra_capacity, length(remainder_non_empty))])
        end
    end

    return keep
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
        # Empty routes are reusable, active routes are unique by key.
        filtered = [e for e in compatible_routes if Int(e.plane.flightLegs) == 0 || !(single_route_key(e.plane) in selected_keys)]
        if isempty(filtered)
            cand_planes[pos] = planeSolution(Int32(0), Int32[base_port], Int32[])
            continue
        end

        # Randomly pick from the top 10 (or fewer if not available)
        # top_k = min(100, length(filtered))
        choice = filtered[rand(1:length(filtered))]

        # Apply chosen route and mark its key as used
        key = single_route_key(choice.plane)
        if Int(choice.plane.flightLegs) > 0
            push!(selected_keys, key)
        end
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


function HeuristicSA(maxTurnaround::Int64, MaxTime::Int32, data, rt, top_c)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0
    iterations_since_clear = 0
    clear_interval = 2000

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    candiateroutes = Candidate_Route(data)
    single_route_pool = SingleRoutePoolEntry[]


    max_size = Int(round(length(data.N)*length(data.V)))*2

    
    while elapsed <= Float64(MaxTime)
        nr = 5
        idx = rand(1:nr)
        T = 100
        Best_sols = generate_best_initial_solutions(data, rt, candiateroutes; n_runs = 10, top_k = nr, top_c, maxLegs=6, maxTurnaround)

        for (i, best_sol_candidate) in enumerate(Best_sols)
            check_route_endpoints(best_sol_candidate.evtols, data, "Initial solution $(i)")
        end

        temp_obj = Best_sols[idx].fitness
        temp_sol = deepcopy(Best_sols[idx].evtols)

        for best_sol_candidate in Best_sols
            single_route_pool = collect_feasible_single_plane_routes(best_sol_candidate.evtols, data, rt; pool=single_route_pool, max_size=max_size, debug=true)
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
                check_route_endpoints(pool_sol, data, "Pool candidate")
                pool_obj = obj(pool_sol, data, rt)
                if pool_obj > temp_obj || rand() < 0.1
                    temp_sol = pool_sol
                    temp_obj = pool_obj
                    single_route_pool = collect_feasible_single_plane_routes(temp_sol, data, rt; pool=single_route_pool, max_size=max_size)
                    check_route_endpoints(temp_sol, data, "Accepted pool candidate")
                end
            end
        end

        improvement = true
        while improvement && elapsed <= Float64(MaxTime)
            improvement = false
            method_used = 0

            cand_obj = temp_obj
            cand_sol = deepcopy(temp_sol)
            n_perm = 3

            for _ in 1:n_perm
                move_choice = rand(1:3)
                if move_choice == 1
                    cand_obj, cand_sol = DestructSA(cand_sol, maxTurnaround, cand_obj, data, rt)
                    next_method = 0
                elseif move_choice == 2
                    cand_obj, cand_sol = ConstructSA(cand_sol, maxTurnaround, cand_obj, data, rt)
                    next_method = 1
                else
                    cand_obj, cand_sol = SegmentExchange(cand_sol, maxTurnaround, cand_obj, data, rt)
                    next_method = 2
                end
                method_used = next_method
                check_route_endpoints(cand_sol, data, "Candidate after $(("DestructSA", "ConstructSA", "SegmentExchange")[next_method + 1])")
            end

            # out = 0
            
            # if cand_obj <= -900000
                    
            #     out = Repair(cand_sol, data, rt, cand_obj)
            #     if isnothing(out)
            #         continue
            #     end
            #     cand_obj = obj(cand_sol, data, rt)

            # end 

            delta = cand_obj - temp_obj
        
            if cand_obj > temp_obj || (rand() < exp(max(min(delta/T, 700.0), -700.0)))
                temp_obj = cand_obj
                temp_sol = deepcopy(cand_sol)
                check_route_endpoints(temp_sol, data, "Accepted candidate")
                improvement = true
                single_route_pool = collect_feasible_single_plane_routes(temp_sol, data, rt; pool=single_route_pool, max_size=max_size)

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                    println("New best Obj $(best_obj)")
                    println("Method used: $(("Destructor", "Constructor", "SegmentExchange")[method_used + 1])")
                end
            end

            elapsed = (time_ns() - start_ns) / 1e9
            T = max(T *0.9, 0.0001)
        end

        iterations += 1
        iterations_since_clear += 1
        if iterations_since_clear >= clear_interval
            sort!(single_route_pool, by = e -> (-e.diversity_score, -e.score))
            keep_n = max(1, cld(length(single_route_pool), 2))
            resize!(single_route_pool, keep_n)
            iterations_since_clear = 0
            println("INFO: Pruned single_route_pool to top $(keep_n) entries at iteration $(iterations)")
        end
        if iterations % 100 == 0
            println("iteration: $(iterations), Time in seconds: $(round(elapsed))")
            println("Obj In this iteration $(temp_obj)")
        end
        elapsed = (time_ns() - start_ns) / 1e9
    end

    println("\n========== FINAL POOL SOLUTIONS ==========")
    println("Total pool entries: $(length(single_route_pool))")
    println("Legend: [idx] flightLegs=N, route=(...), diversity_score=D, fitness_score=F, passengers=[...]")
    println()
    
    # Group by base vertiport for better readability
    entries_by_base = Dict{Int, Vector{Tuple{Int, SingleRoutePoolEntry}}}()
    for (idx, entry) in enumerate(single_route_pool)
        if Int(entry.plane.flightLegs) > 0
            base = Int(entry.plane.route[1])
            if !haskey(entries_by_base, base)
                entries_by_base[base] = []
            end
            push!(entries_by_base[base], (idx, entry))
        else
            if !haskey(entries_by_base, -1)
                entries_by_base[-1] = []
            end
            push!(entries_by_base[-1], (idx, entry))
        end
    end
    
    for (base, entries) in sort(collect(entries_by_base), by=x->x[1])
        if base == -1
            println("--- Empty Routes (Inactive Planes) ---")
        else
            println("--- Base Vertiport $base ---")
        end
        for (idx, entry) in entries
            key = single_route_key(entry.plane)
            passenger_groups = [ass.group for ass in entry.assignments]
            groups_str = isempty(passenger_groups) ? "none" : join(passenger_groups, ",")
            diversity_pct = round(entry.diversity_score * 100; digits=1)
            println("  [$idx] flightLegs=$(key[1]), route=$(key[2]), div=$(diversity_pct)%, fitness=$(round(entry.score; digits=1)), passengers=[$(groups_str)]")
        end
        println()
    end
    
    println("=========================================\n")

    return best_obj, best_sol, iterations
end




