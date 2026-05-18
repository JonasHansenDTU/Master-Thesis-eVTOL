function DestructSA(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    N = data.N

    n = rand(N)
    m = Int(planes.planes[n].flightLegs)

    if m == 0
        return init_obj, planes
    end

    if m == 1
        temp_sol = deepcopy(planes)
        base = temp_sol.planes[n].route[1]
        resize!(temp_sol.planes[n].route, 1)
        temp_sol.planes[n].route[1] = base
        empty!(temp_sol.planes[n].turnaroundTime)
        temp_sol.planes[n].flightLegs = Int32(0)

        new_obj = obj(temp_sol, data, rt)
        return new_obj, temp_sol
    end

    # Do not delete first node (route[1]); last route node is not in this range anyway
    
    idx = rand(2:m)
    temp_sol = deepcopy(planes)

    Destructor(temp_sol.planes[n], [idx], data)

    if has_consecutive_duplicates(temp_sol.planes[n].route)
        return init_obj, planes
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
    
    

    # for n in N
    #     if best_sol.planes[n].route[1] != best_sol.planes[n].route[end]
    #         println("HEY!!!")
    #     end
    # end

    return new_obj, temp_sol
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
            # Fallback: build 2-leg route returning to base
            temp_sol = deepcopy(planes)
            te32 = Int32(data.te)
            insert!(temp_sol.planes[n].route, 2, Int32(base))
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
