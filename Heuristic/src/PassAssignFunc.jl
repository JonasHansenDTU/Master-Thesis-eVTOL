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
                    leg2.to == dp[a] &&
                    leg2.remaining_capacity >= q[a]

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

function export_solution_snapshots(evtols::allPlaneSolution, scheduled::Vector{ScheduledLeg}, assignments::Vector{PassengerAssignment}, data; out_csv::String = joinpath(@__DIR__, "..", "solution_snapshots.csv"))
    V = data.V
    A = data.A
    N = length(evtols.planes)
    lat = data.lat
    lon = data.lon
    rt = data.rt
    op = data.op
    dp = data.dp
    q = data.q
    ET = Int(data.ET)

    # Build a map of which passenger groups are assigned to which legs
    passengers_on_leg = Dict{Tuple{Int,Int}, Vector{Int}}()  # (plane, leg_index) -> [group_ids]
    for ass in assignments
        for leg_idx in ass.legs
            key = (ass.plane, leg_idx)
            if !haskey(passengers_on_leg, key)
                passengers_on_leg[key] = Int[]
            end
            push!(passengers_on_leg[key], ass.group)
        end
    end

    # Build a map of which groups are served overall by each eVTOL
    served_groups_by_evtol = Dict{Int, Vector{Int}}()
    for n in 1:N
        groups = Int[]
        for (plane, leg_idx) in keys(passengers_on_leg)
            if plane == n
                append!(groups, passengers_on_leg[(plane, leg_idx)])
            end
        end
        served_groups_by_evtol[n] = unique(groups)
    end

    rows = NamedTuple[]

    for n in 1:N, t in 0:ET
        # Determine state and location of eVTOL n at time t
        state = "inactive"
        node_from = missing
        node_to = missing
        op_index = missing
        is_p = 0.0
        is_o = 0.0

        # Find which leg(s) this eVTOL is on at time t
        for leg in scheduled
            if leg.plane != n
                continue
            end

            if leg.dep <= t < leg.arr
                # Flying
                state = "flying"
                node_from = leg.from
                node_to = leg.to
                op_index = leg.leg_index
                is_o = 1.0
                break
            elseif leg.arr <= t && (leg.leg_index == length([l for l in scheduled if l.plane == n]) || t < [l for l in scheduled if l.plane == n && l.leg_index > leg.leg_index][1].dep)
                # Parked at destination
                state = "parked"
                node_from = leg.to
                node_to = leg.to
                is_p = 1.0
                break
            end
        end

        # Handle parked at base for first time step before first flight
        if state == "inactive"
            first_leg = findfirst(l -> l.plane == n, scheduled)
            if first_leg !== nothing
                if t < scheduled[first_leg].dep
                    state = "parked"
                    # Find base vertiport for this eVTOL
                    base_vp = data.bv[n]
                    node_from = base_vp
                    node_to = base_vp
                    is_p = 1.0
                else
                    state = "inactive"
                end
            else
                state = "parked"
                base_vp = data.bv[n]
                node_from = base_vp
                node_to = base_vp
                is_p = 1.0
            end
        end

        # Get coordinates
        x_from = ismissing(node_from) ? missing : lon[node_from]
        y_from = ismissing(node_from) ? missing : lat[node_from]
        x_to = ismissing(node_to) ? missing : lon[node_to]
        y_to = ismissing(node_to) ? missing : lat[node_to]

        x = (ismissing(x_from) || ismissing(x_to)) ? missing : (x_from + x_to) / 2
        y = (ismissing(y_from) || ismissing(y_to)) ? missing : (y_from + y_to) / 2

        # Compute battery level for time t using physics-based calculation
        # Start with full battery and simulate charging/discharging
        battery_level = Float32(data.bmax)
        current_time = 0
        
        for leg in sort([l for l in scheduled if l.plane == n], by=l->l.dep)
            # Charge during parked time before this leg
            parked_time = leg.dep - current_time
            battery_level = min(battery_level + Float32(parked_time) * Float32(data.ec), Float32(data.bmax))
            
            # During flight of this leg
            if leg.dep <= t < leg.arr
                # Interpolate position during flight based on discharge rate
                time_in_flight = t - leg.dep
                flight_distance = Float32(data.dist[(leg.from, leg.to)]) 
                discharge_rate = Float32(data.battery_per_km) * flight_distance / Float32(leg.arr - leg.dep)
                battery_level = battery_level - discharge_rate * Float32(time_in_flight)
                break
            elseif t >= leg.arr
                # After flight of this leg, apply full discharge
                flight_distance = Float32(data.dist[(leg.from, leg.to)])
                battery_level = battery_level - Float32(data.battery_per_km) * flight_distance
                current_time = leg.arr
            else
                break
            end
        end
        
        # Clamp battery to valid range
        battery_level = max(battery_level, Float32(data.bmin))

        # Get passenger info
        onboard_groups = Int[]
        if state == "flying" && op_index !== missing
            onboard_groups = get(passengers_on_leg, (n, op_index), Int[])
        end
        onboard_passenger_count = sum(q[a] for a in onboard_groups; init=0)
        onboard_group_sizes = isempty(onboard_groups) ? "" : join(["$(a):$(q[a])" for a in onboard_groups], ";")

        served_groups_evtol = served_groups_by_evtol[n]

        push!(rows, (
            time = t,
            evtol_id = n,
            state = state,
            node_from = node_from,
            node_to = node_to,
            op = op_index,
            is_p = is_p,
            is_o = is_o,
            battery_level = battery_level,
            battery_after_op = op_index,
            onboard_passenger_count = onboard_passenger_count,
            onboard_groups = isempty(onboard_groups) ? "" : join(onboard_groups, ";"),
            onboard_group_sizes = onboard_group_sizes,
            served_groups_evtol = isempty(served_groups_evtol) ? "" : join(served_groups_evtol, ";"),
            x = x,
            y = y,
            x_from = x_from,
            y_from = y_from,
            x_to = x_to,
            y_to = y_to
        ))
    end

    # Append infrastructure nodes as vertiport markers
    for j in V
        push!(rows, (
            time = -1,
            evtol_id = -1,
            state = "vertiport",
            node_from = j,
            node_to = j,
            op = -1,
            is_p = 0.0,
            is_o = 0.0,
            battery_level = NaN,
            battery_after_op = -1,
            onboard_passenger_count = 0,
            onboard_groups = "",
            onboard_group_sizes = "",
            served_groups_evtol = "",
            x = lon[j],
            y = lat[j],
            x_from = lon[j],
            y_from = lat[j],
            x_to = lon[j],
            y_to = lat[j]
        ))
    end

    snapshots = DataFrame(rows)
    CSV.write(out_csv, snapshots)
    println("Snapshot export written: ", out_csv, " (rows=", nrow(snapshots), ")")
    return snapshots
end
