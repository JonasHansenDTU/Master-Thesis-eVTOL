mutable struct PassengerAssignment

    group::Int
    plane::Int
    legs::Vector{Int}
end

mutable struct ScheduledLeg
    plane::Int
    leg_index::Int
    from::Int
    to::Int
    dep::Int
    arr::Int
    remaining_capacity::Int
end

function build_scheduled_legs(evtols::allPlaneSolution, rt::Matrix{Int}, cap_u::Int)
    scheduled = ScheduledLeg[]

    for (pidx, plane) in enumerate(evtols.planes)
        current_time = 0

        for i in 1:plane.flightLegs
            from = Int(plane.route[i])
            to = Int(plane.route[i + 1])

            dep = current_time + Int(plane.turnaroundTime[i])
            arr = dep + rt[from, to]

            push!(scheduled, ScheduledLeg(
                pidx,
                i,
                from,
                to,
                dep,
                arr,
                cap_u
            ))

            current_time = arr
        end
    end

    return scheduled
end

function print_schedule_pretty(scheduled::Vector{ScheduledLeg}; sort_by_plane::Bool=true)
    if isempty(scheduled)
        println("Schedule is empty.")
        return
    end

    legs = sort_by_plane ? sort(scheduled, by = l -> (l.plane, l.leg_index)) : scheduled

    println("Schedule")
    println("========")

    current_plane = -1
    for leg in legs
        if leg.plane != current_plane
            current_plane = leg.plane
            println()
            println("Plane $(current_plane)")
            println("-" ^ 72)
            @printf("%-6s %-6s %-6s %-6s %-8s %-8s %-8s\n",
                    "Leg", "From", "To", "Dep", "Arr", "Dur", "CapLeft")
            println("-" ^ 72)
        end

        dur = leg.arr - leg.dep
        @printf("%-6d %-6d %-6d %-6d %-8d %-8d %-8d\n",
                leg.leg_index, leg.from, leg.to, leg.dep, leg.arr, dur, leg.remaining_capacity)
    end

    println()
    println("Total legs: $(length(legs))")
end

mutable struct PassengerAssignment

    group::Int
    plane::Int
    legs::Vector{Int}
end

function find_direct_leg!(scheduled::Vector{ScheduledLeg}, a::Int, op, dp, dt, q, w)
    candidates = ScheduledLeg[]

    earliest = Int(round(dt[a]))
    latest = earliest + Int(round(w))

    for leg in scheduled
        if leg.from == op[a] &&
           leg.to == dp[a] &&
           earliest <= leg.dep <= latest &&
           leg.remaining_capacity >= q[a]
            push!(candidates, leg)
        end
    end

    if isempty(candidates)
        return nothing
    end

    best = candidates[argmin([leg.arr for leg in candidates])]
    best.remaining_capacity -= q[a]

    return PassengerAssignment(a, best.plane, [best.leg_index])
end

function find_direct_leg!V2(scheduled::Vector{ScheduledLeg}, a::Int, op, dp, dt, q, w, so)
    candidates = ScheduledLeg[]

    earliest = Int(round(dt[a]))
    latest = earliest + Int(round(w))

    for leg in scheduled
        if leg.from == op[a] &&
           leg.to == dp[a] &&
           earliest <= leg.dep <= latest &&
           leg.remaining_capacity >= q[a]
            push!(candidates, leg)
        end
    end

    if isempty(candidates)
        return nothing
    end

    best = candidates[argmin([leg.arr for leg in candidates])]
    if so[a] == 1
        best.remaining_capacity -= q[a]
    else 
        best.remaining_capacity = 0 #Makes sure that the direct only passengers ride alone
    end

    return PassengerAssignment(a, best.plane, [best.leg_index])
end

function find_one_stop_assignment!(scheduled::Vector{ScheduledLeg},
                                   a::Int, op, dp, dt, q, w)
    best_pair = nothing
    best_arrival = typemax(Int)

    earliest = Int(round(dt[a]))
    latest = earliest + Int(round(w))

    for leg1 in scheduled
        if leg1.from != op[a]
            continue
        end
        if !(earliest <= leg1.dep <= latest)
            continue
        end
        if leg1.remaining_capacity < q[a]
            continue
        end

        for leg2 in scheduled  #Kig kun på det næste leg fra leg1
            if leg2.plane != leg1.plane
                continue
            end
            if leg2.leg_index <= leg1.leg_index
                continue
            end
            if leg2.from != leg1.to
                continue
            end
            if leg2.to != dp[a]
                continue
            end
            if leg2.remaining_capacity < q[a]
                continue
            end
            if leg2.dep < leg1.arr
                continue
            end

            if leg2.arr < best_arrival
                best_arrival = leg2.arr
                best_pair = (leg1, leg2)
            end
        end
    end

    if best_pair === nothing
        return nothing
    end

    leg1, leg2 = best_pair
    leg1.remaining_capacity -= q[a]
    leg2.remaining_capacity -= q[a]

    return PassengerAssignment(a, leg1.plane, [leg1.leg_index, leg2.leg_index])
end

function find_one_stop_assignment!V2(scheduled::Vector{ScheduledLeg},
                                   a::Int, op, dp, dt, q, w)
    best_pair = nothing
    best_arrival = typemax(Int)

    earliest = Int(round(dt[a]))
    latest = earliest + Int(round(w))

    for (idx,leg1) in enumerate(scheduled)

        if idx < size(scheduled)[1]
            if leg1.from == op[a] &&
                (earliest <= leg1.dep <= latest)&&
                leg1.remaining_capacity >= q[a]

                leg2 = scheduled[idx+1]
                if leg2.plane == leg1.plane &&
                    leg2.to == dp[a]

                    if leg2.arr < best_arrival
                        best_arrival = leg2.arr
                        best_pair = (leg1, leg2)
                    end
                end
            end
        end
    end

    if best_pair === nothing
        return nothing
    end

    leg1, leg2 = best_pair
    leg1.remaining_capacity -= q[a]
    leg2.remaining_capacity -= q[a]

    return PassengerAssignment(a, leg1.plane, [leg1.leg_index, leg2.leg_index])
end

function assign_passengers(evtols::allPlaneSolution, data, rt::Matrix{Int})
    A = data.A
    op = data.op
    dp = data.dp
    dt = data.dt
    q = data.q
    so = data.so
    w = data.w
    cap_u = Int(round(data.cap_u))

    scheduled = build_scheduled_legs(evtols, rt, cap_u)
    assignments = PassengerAssignment[]
    assigned_groups = Set{Int}()

    # Step 1: non-stopover passengers -> direct only
    direct_only = sort([a for a in A if so[a] == 0], by = a -> (dt[a], -q[a]))
    for a in direct_only
        ass = find_direct_leg!(scheduled, a, op, dp, dt, q, w)
        if ass !== nothing
            push!(assignments, ass)
            push!(assigned_groups, a)
        end
    end

    # Step 2: stopover-allowed passengers -> direct first
    stopover_ok = sort([a for a in A if so[a] == 1], by = a -> (dt[a], -q[a]))
    for a in stopover_ok
        if a in assigned_groups #Not needed?
            continue
        end
        ass = find_direct_leg!(scheduled, a, op, dp, dt, q, w)
        if ass !== nothing
            push!(assignments, ass)
            push!(assigned_groups, a)
        end
    end

    # Step 3: remaining stopover-allowed passengers -> one-stop
    for a in stopover_ok
        if a in assigned_groups
            continue
        end
        ass = find_one_stop_assignment!(scheduled, a, op, dp, dt, q, w)
        if ass !== nothing
            push!(assignments, ass)
            push!(assigned_groups, a)
        end
    end

    return assignments, scheduled
end

function assign_passengersV2(evtols::allPlaneSolution, data, rt::Matrix{Int})
    A = data.A
    op = data.op
    dp = data.dp
    dt = data.dt
    q = data.q
    so = data.so
    w = data.w
    fd = data.fd
    fs = data.fs
    cap_u = Int(round(data.cap_u))

    scheduled = build_scheduled_legs(evtols, rt, cap_u)
    assignments = PassengerAssignment[]
    assigned_groups = Set{Int}()


    price_by_group = Dict(a => fd[(op[a], dp[a])] * (so[a] == 1 ? 0.75 : 1.0) for a in A)
    Price_sort = sort(A, by = a -> price_by_group[a], rev = true)  #descending


    for a in Price_sort
        if so[a] == 0
            ass = find_direct_leg!V2(scheduled, a, op, dp, dt, q, w, so)
            if ass !== nothing
                push!(assignments, ass)
                push!(assigned_groups, a)
            end
        else
            ass = find_direct_leg!V2(scheduled, a, op, dp, dt, q, w, so)
            if ass == nothing 
                ass = find_one_stop_assignment!V2(scheduled, a, op, dp, dt, q, w)
                if ass !== nothing
                    push!(assignments, ass)
                    push!(assigned_groups, a)
                end
            else
                push!(assignments, ass)
                push!(assigned_groups, a)
            end
        end
    end

    return assignments, scheduled
end

function print_assignments(assignments::Vector{PassengerAssignment}, data)
    println("Passenger assignments:")
    println("----------------------")

    for ass in assignments
        println(
            "Group ", ass.group,
            " | plane ", ass.plane,
            " | legs ", ass.legs,
            " | ", data.op[ass.group], " -> ", data.dp[ass.group],
            " | q=", data.q[ass.group],
            " | dt=", data.dt[ass.group]
        )
    end
end
