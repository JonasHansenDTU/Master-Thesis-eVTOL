function solution_kpi_context(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    rt_int = Int.(rt)
    assignments, scheduled = assign_passengersV2(evtols, data, rt_int)

    passenger_map = Dict{Tuple{Int,Int}, Vector{Int}}()
    for ass in assignments
        for leg_idx in ass.legs
            key = (ass.plane, leg_idx)
            if !haskey(passenger_map, key)
                passenger_map[key] = Int[]
            end
            push!(passenger_map[key], ass.group)
        end
    end

    scheduled_lookup = Dict((leg.plane, leg.leg_index) => leg for leg in scheduled)

    return assignments, scheduled, passenger_map, scheduled_lookup
end

function aircraft_utilization(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    # Total available capacity over the time horizon (Total Planes * Capacity * Total Time)
    # If data.ET is just time, we scale it by the individual plane capacity
    total_planes = max(length(evtols.planes), 1)
    active_planes = 0

    for i in 1:total_planes
        if evtols.planes[i].flightLegs > 1
            active_planes += 1
        end
    end



    total_capacity_available = active_planes * Float64(data.cap_u) * Float64(data.ET)
    
    total_capacity_occupied = 0.0
    minimum_turnaround = Float64(data.te)
    assignments, scheduled, passenger_map, scheduled_lookup = solution_kpi_context(evtols, data, rt)

    for (plane_idx, plane) in enumerate(evtols.planes)
        flight_seat_hours = 0.0
        turnaround_seat_hours = 0.0

        for leg_idx in 1:plane.flightLegs
            from = Int(plane.route[leg_idx])
            to = Int(plane.route[leg_idx + 1])            
            leg = scheduled_lookup[(plane_idx, leg_idx)]

            # 1. Calculate actual duration of the flight leg
            leg_duration = Float64(rt[from, to])
            
            # 2. Calculate passengers on board (Total capacity minus what is left)
            passengers_on_board = Float64(data.cap_u - leg.remaining_capacity)

            # 3. Scale the flight duration by actual passenger load factor
            # e.g., 0.5 hours with 3 passengers = 1.5 passenger-hours
            flight_seat_hours += leg_duration * passengers_on_board
            
            # 4. Turnaround utilizes the full asset capacity for operational overhead
            turnaround_seat_hours += minimum_turnaround * Float64(data.cap_u)
        end

        total_capacity_occupied += flight_seat_hours + turnaround_seat_hours
    end

    return total_capacity_available == 0 ? 0.0 : total_capacity_occupied / total_capacity_available
end

function deadhead_time(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    _, scheduled, passenger_map, _ = solution_kpi_context(evtols, data, rt)

    deadhead = 0.0
    for leg in scheduled
        if leg.remaining_capacity == data.cap_u
            deadhead += Float64(leg.arr - leg.dep)
        end
    end

    return deadhead
end

function passenger_demand_served(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    assignments, _, _, _ = solution_kpi_context(evtols, data, rt)

    served_demand = sum(1 for ass in assignments; init=0)

    return served_demand
end

function average_passenger_waiting_time(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    assignments, scheduled, _, scheduled_lookup = solution_kpi_context(evtols, data, rt)

    if isempty(assignments)
        return 0.0
    end

    total_wait = 0.0
    served_groups = 0

    for ass in assignments
        leg_deps = Float64[]
        for leg_idx in ass.legs
            leg = scheduled_lookup[(ass.plane, leg_idx)]
            push!(leg_deps, Float64(leg.dep))
        end

        first_departure = minimum(leg_deps)
        total_wait += max(first_departure - Float64(data.dt[ass.group]), 0.0)
        served_groups += 1
    end

    return total_wait / served_groups
end

function number_of_stops_per_hour(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    total_stops = sum(Float64(plane.flightLegs) for plane in evtols.planes)
    total_available_hours = max(length(evtols.planes), 1) * Float64(data.ET) / 60.0 # THIS IS IMPORTANT IF UNITTIME IS CHANGED!!!!

    return total_available_hours == 0 ? 0.0 : total_stops / total_available_hours
end

function aircraft_idle_time(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    total_idle = 0.0
    minimum_turnaround = Float64(data.te)

    for plane in evtols.planes
        for leg_idx in 1:plane.flightLegs
            extra_turnaround = Float64(plane.turnaroundTime[leg_idx]) - minimum_turnaround
            total_idle += max(extra_turnaround, 0.0)
        end
    end

    return isempty(evtols.planes) ? 0.0 : total_idle / length(evtols.planes)
end

function solution_kpis(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    return (
        aircraft_utilization = aircraft_utilization(evtols, data, rt),
        deadhead_time = deadhead_time(evtols, data, rt),
        passenger_demand_served = passenger_demand_served(evtols, data, rt),
        average_passenger_waiting_time = average_passenger_waiting_time(evtols, data, rt),
        # number_of_stops_per_hour = number_of_stops_per_hour(evtols, data, rt),
        # aircraft_idle_time = aircraft_idle_time(evtols, data, rt),
    )
end