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
    total_occupied = 0.0
    total_available = max(length(evtols.planes), 1) * Float64(data.ET)
    minimum_turnaround = Float64(data.te)
    assignments, scheduled, passenger_map, scheduled_lookup = solution_kpi_context(evtols, data, rt)

    for plane in evtols.planes
        flight_time = 0.0
        occupied_turnaround = 0.0

        for leg_idx in 1:plane.flightLegs
            from = Int(plane.route[leg_idx])
            to = Int(plane.route[leg_idx + 1])            
            leg = scheduled_lookup[(plane, leg_idx)]'


            flight_time += Float64(rt[from, to])/(data.cap_u - leg.remaining_capacity)
            occupied_turnaround += minimum_turnaround
        end

        total_occupied += flight_time + occupied_turnaround
    end

    return total_available == 0 ? 0.0 : total_occupied / total_available
end

function deadhead_time(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    _, scheduled, passenger_map, _ = solution_kpi_context(evtols, data, rt)

    deadhead = 0.0
    for leg in scheduled
        key = (leg.plane, leg.leg_index)
        if !haskey(passenger_map, key) || isempty(passenger_map[key])
            deadhead += Float64(leg.arr - leg.dep)
        end
    end

    return deadhead
end

function passenger_demand_served(evtols::allPlaneSolution, data, rt::AbstractMatrix)
    assignments, _, _, _ = solution_kpi_context(evtols, data, rt)

    served_demand = sum(Float64(data.q[ass.group]) for ass in assignments)

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