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

    
    while elapsed <= Float64(MaxTime)
        nr = 5
        idx = rand(1:nr)
        T = 50
        Best_sols = generate_best_initial_solutions(data, rt, candiateroutes; n_runs = 10, top_k = nr, top_c, maxLegs=6, maxTurnaround)

        temp_obj = Best_sols[idx].fitness
        temp_sol = deepcopy(Best_sols[idx].evtols)

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

            out = 0
            
            if cand_obj <= -900000
                    
                out = Repair(cand_sol, data, rt, cand_obj)
                if isnothing(out)
                    continue
                end
                cand_obj = obj(cand_sol, data, rt)

            end 


            if cand_obj > temp_obj || (rand() < exp((temp_obj - cand_obj) / T) && temp_obj > cand_obj)
                temp_obj = cand_obj
                temp_sol = deepcopy(cand_sol)
                improvement = true

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

    return best_obj, best_sol, iterations
end




