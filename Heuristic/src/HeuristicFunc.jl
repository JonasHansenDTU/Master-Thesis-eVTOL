function Heuristic(maxTurnaround::Int64, MaxTime::Int32, data, rt, top_c)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    candiateroutes = Candidate_Route(data)

    while elapsed <= Float64(MaxTime)
        nr = 1
        Best_sols = generate_best_initial_solutions(data, rt, candiateroutes; n_runs = 50, top_k = nr, top_c, maxLegs=6, maxTurnaround)

        temp_obj = Best_sols[nr].fitness
        temp_sol = deepcopy(Best_sols[nr].evtols)

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

            cand_obj = des_obj
            cand_sol = des_sol

            if con_obj > cand_obj
                cand_obj = con_obj
                cand_sol = con_sol
                method_used = 1
            end

            if swap_obj > cand_obj
                cand_obj = swap_obj
                cand_sol = swap_sol
                method_used = 2
            end

            if two_opt_obj > cand_obj
                cand_obj = two_opt_obj
                cand_sol = two_opt_sol
                method_used = 3
            end

            if cand_obj > temp_obj
                temp_obj = cand_obj
                temp_sol = deepcopy(cand_sol)
                improvement = true

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                    println("New best Obj $(best_obj)")
                    println("Method used: $(("Destructor", "Constructor", "Swap", "2-opt")[method_used + 1])")
                end
            end

            elapsed = (time_ns() - start_ns) / 1e9
        end

        iterations += 1
        if iterations % 100 == 0
            println("iteration: $(iterations)")
        end
        elapsed = (time_ns() - start_ns) / 1e9
    end

    return best_obj, best_sol, iterations
end