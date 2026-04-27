function fitnessFunction(evtols::allPlaneSolution, assignments::Vector{PassengerAssignment}, bmax::Float32, bmid::Float32, bmin::Float32,
    dist::Dict{Tuple{Int,Int},Float64}, ec::Float32, battery_per_km::Float32, rt::Matrix{Int}, ET::Int, T::Int, V::Int, cap_v::Dict{}, data)
    
    b_penalty = data.b_penalty 
    A  = data.A
    op = data.op
    dp = data.dp
    fd = data.fd
    fs = data.fs
    c  = data.c
    so = data.so
    p  = data.p
    opening_cost = data.opening_cost

    P = FeasibilityCheck(bmax, bmid, bmin, dist, ec, battery_per_km, evtols, rt, ET, T, V, cap_v, b_penalty)

    fitnessvalue = 0.0

    # 1. Revenue from assigned passenger groups
    for ass in assignments
        a = ass.group
        i = op[a]
        j = dp[a]
        fitnessvalue += fd[(i, j)] * (1 - so[a]) + fs[(i, j)] * so[a]
    end

    # 2. Operating cost of all flown legs
    for plane in evtols.planes
        for k in 1:plane.flightLegs
            from = Int(plane.route[k])
            to   = Int(plane.route[k + 1])
            fitnessvalue -= c[(from, to)]
        end
    end

    # 3. Penalty for unserved passenger groups
    assigned_groups = Set(ass.group for ass in assignments)
    for a in A
        if !(a in assigned_groups)
            fitnessvalue -= p[a]
        end
    end

    # 4. Fixed cost for each used eVTOL
    for plane in evtols.planes
        if plane.flightLegs > 0
            fitnessvalue -= opening_cost
        end
    end

    # 5. Large infeasibility penalty
    fitnessvalue -= sum(P)

    return fitnessvalue
end
