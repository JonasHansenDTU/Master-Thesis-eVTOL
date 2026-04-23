
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

function initial_chromosome_solution(data; maxLegs::Int=5, maxTurnaround::Int=30)
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


    weights = build_vertiport_weights(V, op, dp, A)

    planes = planeSolution[]
    nPlanes = length(N)

    choices = [i for i in 0:maxLegs*nPlanes if i != 1]
    total_flight_legs = rand(choices)
    Set_of_flightlegs = [i for i in 0:total_flight_legs if i != 1 && i != total_flight_legs-1]

    allowed_legs = zeros(Int, nPlanes+1)
    allowed_legs[nPlanes+1] = total_flight_legs

    # println(Set_of_flightlegs)

    for n in 1:nPlanes-1
        allowed_legs[n+1] = (rand(Set_of_flightlegs))
        
        while n > 1 && (allowed_legs[n+1] == allowed_legs[n] -1 
                    || allowed_legs[n+1] == allowed_legs[n] +1) 
            allowed_legs[n] = rand(Set_of_flightlegs) 
            # println("HERE")
        end
    end

    allowed_legs = sort(allowed_legs)


    allowed_legs = shuffle!([allowed_legs[i]-allowed_legs[i-1] for i in 2:nPlanes+1])


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
            from_current = [(a, v) for (a, v) in op if v == current_VP]
            Candidate_Pass = [dt[a] for (a, _) in from_current if current_time+te <= dt[a] <= current_time + maxTurnaround]
            if !isempty(Candidate_Pass)
                chosen_time = rand(Candidate_Pass) - current_time
                push!(turnaroundTime, round(Int32, chosen_time))
            else
                push!(turnaroundTime, Int32(rand(te:maxTurnaround)))
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
    # Opt_evtol2 = []
    # Opt_evtol3 = [5,2,3,5]

    # if planes[1].route == Opt_evtol1 && 
    #     planes[2].flightLegs == 0 &&
    #     planes[3].route == Opt_evtol3

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


    #-----------------------------------------


    return allPlaneSolution(planes)
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


function generate_best_initial_solutions(data, rt; n_runs::Int=1000, top_k::Int=10, maxLegs::Int=5, maxTurnaround::Int=30, print_each::Bool=false)
    results = NamedTuple[]

    for run in 1:n_runs
        evtols_init = initial_chromosome_solution(data; maxLegs=maxLegs, maxTurnaround=maxTurnaround)

        assignments, scheduled = assign_passengersV2(evtols_init, data, rt)

        # model, scheduled = assign_passengers_Solver2(evtols_init, data, rt)
        # assignments = extract_assignments2(model)

        P = FeasibilityCheck(
            Float32(data.bmax),
            Float32(data.bmin),
            data.dist,
            Float32(data.ec),
            Float32(data.battery_per_km),
            evtols_init,
            rt,
            Int(round(data.ET)),
            maximum(Int.(data.T)),
            maximum(data.V),
            Int(round(data.cap_flt)),
            data.cap_v
        )

        fitness = fitnessFunction(
            evtols_init,
            assignments,
            Float32(data.bmax),
            Float32(data.bmin),
            data.dist,
            Float32(data.ec),
            Float32(data.battery_per_km),
            rt,
            Int(round(data.ET)),
            maximum(Int.(data.T)),
            maximum(data.V),
            Int(round(data.cap_flt)),
            data.cap_v,
            data
        )

        push!(results, (
            run = run,
            fitness = fitness,
            evtols = evtols_init,
            assignments = assignments,
            scheduled = scheduled,
            P = P
        ))

        if print_each
            println("Run $run | fitness = $fitness | P = $P")
        end
    end

    # Sort by fitness descending, since higher fitness is better
    sorted_results = sort(results, by = x -> x.fitness, rev = true)

    # Return only the best top_k results
    return sorted_results[1:min(top_k, length(sorted_results))]
end
