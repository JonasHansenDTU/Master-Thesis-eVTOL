using XLSX
using DataFrames
using Random
using JuMP
using Gurobi
using CSV
using MathOptInterface
using Printf
const MOI = MathOptInterface

function normalize_name(x)
    s = lowercase(strip(String(x)))
    s = replace(s, " " => "_")
    s = replace(s, "-" => "_")
    s = replace(s, "/" => "_")
    s = replace(s, "." => "_")
    s = replace(s, "(" => "")
    s = replace(s, ")" => "")
    return Symbol(s)
end

"Read one Excel sheet into a DataFrame with normalized column names."
function read_sheet(path::String, sheet_name::String)
    df = DataFrame(XLSX.readtable(path, sheet_name, infer_eltypes=true))
    rename!(df, Dict(n => normalize_name(n) for n in names(df)))
    return df
end

"Read first available sheet among candidates into a DataFrame."
function read_sheet_any(path::String, candidates::Vector{String})
    for sheet_name in candidates
        try
            return read_sheet(path, sheet_name)
        catch
            # Try next candidate sheet name.
        end
    end

    error("Could not find any of these sheets: $(candidates)")
end

"Try to find a column among several possible names."
function find_col(df::DataFrame, candidates::Vector{Symbol})
    name_map = Dict(Symbol(String(n)) => n for n in names(df))
    for c in candidates
        if haskey(name_map, c)
            return name_map[c]
        end
    end
    error("Could not find any of these columns in sheet: $(candidates). Found columns: $(names(df))")
end

"Parse coordinate string like '(37.5655, 126.8013)' -> (lat, lon)."
function parse_coordinate_string(s)
    ss = strip(String(s))
    ss = replace(ss, "(" => "")
    ss = replace(ss, ")" => "")
    parts = split(ss, ",")
    if length(parts) != 2
        error("Could not parse coordinate string: $s")
    end
    lat = parse(Float64, strip(parts[1]))
    lon = parse(Float64, strip(parts[2]))
    return lat, lon
end

"Haversine distance in km."
function haversine_km(lat1, lon1, lat2, lon2)
    R = 6371.0
    φ1 = deg2rad(lat1)
    λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2)
    λ2 = deg2rad(lon2)
    dφ = φ2 - φ1
    dλ = λ2 - λ1
    a = sin(dφ / 2)^2 + cos(φ1) * cos(φ2) * sin(dλ / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return R * c
end

"""
Build all sets and parameters from Excel + system parameters.
Expected sheets:
    - Infrastructure (3)
    - PassengerGroups (3)
    - PlaneData
"""
function load_data(excel_file::String)

    ###########################################################################
    # Read sheets
    ###########################################################################
    infra = read_sheet(excel_file, "Infrastructure (3)")
    pax   = read_sheet(excel_file, "PassengerGroups (3)")
    plane = read_sheet_any(excel_file, ["PlaneData (2)"])

    ###########################################################################
    # Infrastructure columns
    ###########################################################################
    id_col     = find_col(infra, [:id])
    pads_col   = find_col(infra, [:number_of_parking_pads, :parking_pads, :pads])

    # Coordinates can be either one string column "coordinates"
    # or two numeric columns such as "latitude", "longitude".
    coord_col = if any(Symbol(String(n)) == :coordinates_kbh for n in names(infra))
        find_col(infra, [:coordinates_kbh])
    else
        nothing
    end
    
    lat_col = if any(Symbol(String(n)) == :latitude for n in names(infra))
        find_col(infra, [:latitude])
    elseif any(Symbol(String(n)) == :lat for n in names(infra))
        find_col(infra, [:lat])
    else
        nothing
    end
    
    lon_col = if any(Symbol(String(n)) == :longitude for n in names(infra))
        find_col(infra, [:longitude])
    elseif any(Symbol(String(n)) == :lon for n in names(infra))
        find_col(infra, [:lon])
    else
        nothing
    end

    ###########################################################################
    # Passenger columns
    ###########################################################################
    group_col  = find_col(pax, [:group, :group_id, :id])
    orig_col   = find_col(pax, [:origin])
    dest_col   = find_col(pax, [:destination])
    time_col   = find_col(pax, [:time, :arrival_time])
    q_col      = find_col(pax, [:number_of_passengers, :passengers, :numberofpassengers])
    stop_col   = find_col(pax, [:stopover_allowed, :stopoverallowed, :stopover])

    ###########################################################################
    # Sets
    ###########################################################################
    V = sort(unique(Int.(infra[!, id_col])))

    plane_id_col = find_col(plane, [:plane_id, :id, :evtol_id, :plane])
    base_vp_col  = find_col(plane, [:base_vertiport, :base_vertiport_id, :base, :home_vertiport])

    N = sort(Int.(plane[!, plane_id_col]))
    if isempty(N)
        error("PlaneData sheet is empty. Add at least one row with Plane ID and Base Vertiport.")
    end
    if length(unique(N)) != length(N)
        error("PlaneData contains duplicate Plane ID values. Each Plane ID must be unique.")
    end

    bv = Dict{Int,Int}()  # base vertiport for each eVTOL
    for r in eachrow(plane)
        n = Int(r[plane_id_col])
        b = Int(r[base_vp_col])
        bv[n] = b
    end

    
    bad_bases = sort([b for b in values(bv) if !(b in V)])
    if !isempty(bad_bases)
        error("PlaneData has invalid Base Vertiport values $(bad_bases). Valid vertiports from Infrastructure are $(V).")
    end

    M = 0:6
    M_no0 = 1:maximum(M)
    M_mid = 1:(maximum(M)-1)
    M_no_last = 0:(maximum(M)-1)

    T = 0:120
    T_no0 = 1:maximum(T)

    # Passenger groups
    A = sort(Int.(pax[!, group_col]))

    ###########################################################################
    # System parameters from Table 3
    ###########################################################################
    
    df = XLSX.readtable(excel_file, "Parameters")
    params = Dict{String, Float64}()

    XLSX.openxlsx(excel_file) do xf
        sheet = xf["Parameters"]
        
        for row in XLSX.eachrow(sheet)
            key = row[1]
            val = row[2]
            
            if key !== missing && val !== missing
                key_str = String(key)

                # håndter komma som decimalseparator
                val_str = replace(string(val), "," => ".")
                params[key_str] = parse(Float64, val_str)
            end
        end
    end

    fare_direct_per_km    = params["fare_direct_per_km"]
    fare_stopover_factor  = params["fare_stopover_factor"]
    p_penalty             = params["p_penalty"]
    operating_cost_per_km = params["operating_cost_per_km"]
    battery_per_km        = params["battery_per_km"]
    time_per_km           = params["time_per_km"]
    cap_flt               = params["cap_flt"]
    cap_u                 = params["cap_u"]
    bmax                  = params["bmax"]
    bmin                  = params["bmin"]
    ec                    = params["ec"]
    te                    = params["te"]
    w                     = params["w"]
    ET                    = params["ET"]
    M1                    = 1
    M2a                   = 0
    M2b                   = bmax
    M2c                   = bmax + ec * ET
    M3                    = ET


    ###########################################################################
    # Node coordinates and parking capacities
    ###########################################################################
    lat = Dict{Int,Float64}()
    lon = Dict{Int,Float64}()
    cap_v = Dict{Int,Int}()

    for r in eachrow(infra)
        j = Int(r[id_col])
        cap_v[j] = Int(r[pads_col])

        if coord_col !== nothing
            lat[j], lon[j] = parse_coordinate_string(r[coord_col])
        else
            if lat_col === nothing || lon_col === nothing
                error("Infrastructure sheet must contain either a 'coordinates' column or latitude/longitude columns.")
            end
            lat[j] = Float64(r[lat_col])
            lon[j] = Float64(r[lon_col])
        end
    end

    ###########################################################################
    # Derived arc parameters: distance, fd, fs, c, e, rt
    ###########################################################################
    dist = Dict{Tuple{Int,Int},Float64}()
    fd   = Dict{Tuple{Int,Int},Float64}()
    fs   = Dict{Tuple{Int,Int},Float64}()
    c    = Dict{Tuple{Int,Int},Float64}()
    e    = Dict{Tuple{Int,Int},Float64}()
    rt   = Dict{Tuple{Int,Int},Int}()
  

    for i in V, j in V
        dij = haversine_km(lat[i], lon[i], lat[j], lon[j])
        dist[(i,j)] = dij
        fd[(i,j)] = fare_direct_per_km * dij
        fs[(i,j)] = fare_stopover_factor * fd[(i,j)]
        c[(i,j)]  = operating_cost_per_km * dij
        e[(i,j)]  = battery_per_km * dij
        rt[(i,j)] = Int(ceil(time_per_km * dij))
    end

    ###########################################################################
    # Passenger parameters
    ###########################################################################
    op = Dict{Int,Int}()
    dp = Dict{Int,Int}()
    dt = Dict{Int,Float64}()
    q  = Dict{Int,Int}()
    so = Dict{Int,Int}()
    p  = Dict{Int,Float64}()

    for r in eachrow(pax)
        a = Int(r[group_col])
        op[a] = Int(r[orig_col])
        dp[a] = Int(r[dest_col])
        dt[a] = Float64(r[time_col])
        q[a]  = Int(r[q_col])
        so[a] = Int(r[stop_col])
        p[a]  = p_penalty
    end

    ###########################################################################
    # d[a,i,j] = 1 if passenger group a travels from i to j, else 0
    ###########################################################################
    d = Dict{Tuple{Int,Int,Int},Int}()
    for a in A, i in V, j in V
        d[(a,i,j)] = (i == op[a] && j == dp[a]) ? 1 : 0
    end

    return (
        infra = infra,
        pax = pax,
        plane = plane,
        V = V, A = A, N = collect(N),
        M = collect(M), M_no0 = collect(M_no0), M_mid = collect(M_mid), M_no_last = collect(M_no_last),
        T = collect(T), T_no0 = collect(T_no0),
        bv = bv,
        lat = lat, lon = lon,
        dist = dist, fd = fd, fs = fs, c = c, e = e, rt = rt,
        op = op, dp = dp, dt = dt, q = q, so = so, p = p, d = d,
        cap_v = cap_v, cap_flt = cap_flt, cap_u = cap_u,
        bmax = bmax, bmin = bmin, ec = ec, te = te, w = w, ET = ET, M1 = M1, M2a = M2a, M2b = M2b, M2c = M2c, M3 = M3, battery_per_km = battery_per_km
    )
end

mutable struct planeSolution
    flightLegs::Int32

    route::Array{Int32,1}

    turnaroundTime::Array{Int32,1}
end

mutable struct allPlaneSolution
    planes::Vector{planeSolution}
end

function BatteryNeeded(TravelLength::Float32, battery_per_km::Float32)
    return TravelLength * battery_per_km
end

function BatteryCharged(turnaroundTime::Float16, ec::Float32)
    return turnaroundTime*ec
end

function FeasibleBattery(evtols::allPlaneSolution, bmax::Float32, bmin::Float32, dist::Dict{Tuple{Int,Int},Float64}, ec::Float32, battery_per_km::Float32)
    for evtol in evtols.planes
        BatteryLevel = zeros(Float32, evtol.flightLegs + 1)
        BatteryLevel[1] = bmax

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i + 1]
            TravelLength = Float32(dist[(from, to)])

            BatteryLevel[i + 1] =
                min(BatteryLevel[i] + BatteryCharged(Float16(evtol.turnaroundTime[i]), ec), bmax) -
                BatteryNeeded(TravelLength, battery_per_km)

            if BatteryLevel[i + 1] < bmin
                return false
            end
        end
    end 

    return true
end

function FeasibleCompletionTime(evtols::allPlaneSolution, rt::Matrix{Int}, ET::Int)
    for evtol in evtols.planes
        travelTime = 0

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i+1]

            travelTime += evtol.turnaroundTime[i] + rt[from, to]

            if travelTime > ET
                return false
            end
        end
    end

    return true
end

function FeasibleVertiportCapacity(evtols::allPlaneSolution, rt::Matrix{Int}, T::Int, V::Int, cap_v::Dict{Int64, Int64}, ET::Int)
    parkingTimes = zeros(Int, V, T)

    for evtol in evtols.planes
        travelTime = 0

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i+1]

            startTime = travelTime
            endTime = travelTime + evtol.turnaroundTime[i]

            if endTime > ET 
                endTime = ET
            end

            for t in startTime:endTime
                if t >= 1
                    parkingTimes[from, t] += 1
                    if parkingTimes[from, t] > cap_v[from]
                        return false
                    end
                end
            end

            travelTime += evtol.turnaroundTime[i] + rt[from, to]
        end
    end

    return true
end

function FeasibleCorridor(evtols::allPlaneSolution, rt::Matrix{Int}, T::Int, V::Int, cap_flt::Int, ET::Int)
    destinationTimes = zeros(Int, V, V, T)

    for evtol in evtols.planes
        travelTime = 0

        for i in 1:evtol.flightLegs
            from = evtol.route[i]
            to = evtol.route[i+1]

            startTime = travelTime + evtol.turnaroundTime[i]
            endTime = travelTime + evtol.turnaroundTime[i] + rt[from, to]

            if endTime > ET 
                endTime = ET
            end

            for t in startTime:endTime
                if t >= 1
                    destinationTimes[from, to, t] += 1
                    if destinationTimes[from, to, t] > cap_flt 
                        return false
                    end
                end
            end

            travelTime += evtol.turnaroundTime[i] + rt[from, to]
        end
    end

    return true
end

function FeasibilityCheck(bmax::Float32, bmin::Float32,
    dist::Dict{Tuple{Int,Int},Float64},
    ec::Float32, battery_per_km::Float32,
    evtols::allPlaneSolution,
    rt::Matrix{Int}, ET::Int, T::Int, V::Int,
    cap_flt::Int, cap_v::Dict{})

    P = zeros(Int32, 4)

    if FeasibleBattery(evtols, bmax, bmin, dist, ec, battery_per_km) == false
        P[1] = 1
    end

    if FeasibleCompletionTime(evtols, rt, ET) == false
        P[2] = 1
    end

    if FeasibleVertiportCapacity(evtols, rt, T, V, cap_v, ET) == false
        P[3] = 1
    end

    if FeasibleCorridor(evtols, rt, T, V, cap_flt, ET) == false
        P[4] = 1
    end

    return P
end


### ------- TESTING AREA !!! --------------###
# Test if the current version of FeasibleBattery battery works 
# dist = Dict{Tuple{Int,Int},Float64}()
# dist[(1,5)] = 10.0
# dist[(5,3)] = 8.0

# bmax = Float32(100.0)
# bmin = Float32(20.0)
# ec = Float32(1.0)
# battery_per_km = Float32(5.0)

# # true example
# evtol_true = planeSolution(
#     Int32(2),
#     Int32[1, 5, 3],
#     Int32[10, 10]
# )

# # false example
# evtol_false = planeSolution(
#     Int32(2),
#     Int32[1, 5, 3],
#     Int32[0, 0]
# )

# println("True example (should be true):")
# println(FeasibleBattery(evtol_true, bmax, bmin, dist, ec, battery_per_km))

# println("False example (should be false):")
# println(FeasibleBattery(evtol_false, bmax, bmin, dist, ec, battery_per_km))

#Test if this current version of completion time works
# V = 3
# rt = [
#     0 2 1;
#     2 0 3;
#     1 3 0
# ]
# ET = 10

# # false example
# p1 = planeSolution(Int32(2), Int32[1, 2, 3], Int32[5, 5])
# evtols1 = allPlaneSolution([p1])

# println("False example (should be false):")
# println(FeasibleCompletionTime(evtols1, rt, ET))   # false

# # true example
# p2 = planeSolution(Int32(1), Int32[1, 2], Int32[5])
# evtols2 = allPlaneSolution([p2])

# println("True example (should be true):")
# println(FeasibleCompletionTime(evtols2, rt, ET))   # true


#Test if the current version of vertiport capacity works
# V = 3
# T = 20

# rt = [
#     0 2 1;
#     2 0 3;
#     1 3 0
# ]

# cap_v = [1, 2, 2]

# p1 = planeSolution(Int32(1), Int32[1, 2], Int32[5])
# p2 = planeSolution(Int32(1), Int32[1, 3], Int32[5])

# evtols_false = allPlaneSolution([p1, p2])

# println("False example (should be false):")
# println(FeasibleVertiportCapacity(evtols_false, rt, T, V, cap_v))

# p3 = planeSolution(Int32(1), Int32[1, 2], Int32[5])
# p4 = planeSolution(Int32(1), Int32[1, 3], Int32[0])

# evtols_true = allPlaneSolution([p3, p4])

# println("True example (should be true):")
# println(FeasibleVertiportCapacity(evtols_true, rt, T, V, cap_v))


#Test if the current version of corridor works
# V = 3
# T = 20
# cap = 1

# rt = [
#     0 2 1;
#     2 0 3;
#     1 3 0
# ]

# p1 = planeSolution(Int32(1), Int32[1, 2], Int32[1])
# p2 = planeSolution(Int32(1), Int32[1, 2], Int32[1])

# evtols_false = allPlaneSolution([p1, p2])

# println("False example (should be false):")
# println(FeasibleCorridor(evtols_false, rt, T, V, cap))

# p3 = planeSolution(Int32(1), Int32[1, 2], Int32[1])
# p4 = planeSolution(Int32(1), Int32[1, 2], Int32[5])

# evtols_true = allPlaneSolution([p3, p4])

# println("True example (should be true):")
# println(FeasibleCorridor(evtols_true, rt, T, V, cap))

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
            if ass === nothing 
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

function assign_passengers_Solver2(evtols::allPlaneSolution, data, rt::Matrix{Int})
    A  = data.A
    V  = data.V
    N  = 1:length(evtols.planes)
    M  = data.M

    op = data.op
    dp = data.dp
    dt = data.dt
    q  = data.q
    so = data.so
    w  = Int(round(data.w))

    fd = data.fd
    fs = data.fs
    p  = data.p

    cap_u = Int(round(data.cap_u))

    scheduled = build_scheduled_legs(evtols, rt, cap_u)

    ###########################################################################
    # Fixed schedule lookup
    ###########################################################################
    leg_of = Dict{Tuple{Int,Int}, ScheduledLeg}()
    for leg in scheduled
        leg_of[(leg.leg_index, leg.plane)] = leg
    end

    ###########################################################################
    # Feasible assignment keys
    ###########################################################################
    direct_keys = Tuple{Int,Int,Int}[]
    stop_keys = Tuple{Int,Int,Int}[]

    for a in A
        earliest = Int(round(dt[a]))
        latest = earliest + w

        # Direct assignments
        for n in N
            for m in 1:maximum(M)
                if haskey(leg_of, (m, n))
                    leg = leg_of[(m, n)]
                    if leg.from == op[a] &&
                       leg.to == dp[a] &&
                       earliest <= leg.dep <= latest
                        push!(direct_keys, (a, m, n))
                    end
                end
            end
        end

        # One-stop assignments on consecutive legs
        if so[a] == 1
            for n in N
                for m in 1:(maximum(M) - 1)
                    if haskey(leg_of, (m, n)) && haskey(leg_of, (m + 1, n))
                        leg1 = leg_of[(m, n)]
                        leg2 = leg_of[(m + 1, n)]

                        if leg1.from == op[a] &&
                           leg2.from == leg1.to &&
                           leg2.to == dp[a] &&
                           earliest <= leg1.dep <= latest &&
                           leg2.dep >= leg1.arr
                            push!(stop_keys, (a, m, n))
                        end
                    end
                end
            end
        end
    end

    ###########################################################################
    # Solver
    ###########################################################################
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", 0)

    @variable(model, y_dir[direct_keys], Bin)
    @variable(model, y_stop[stop_keys], Bin)

    @objective(
        model, Max,
        sum(
            (fd[(op[a], dp[a])] * (1 - so[a]) + fs[(op[a], dp[a])] * so[a] + p[a]) * y_dir[(a, m, n)]
            for (a, m, n) in direct_keys
        )
        +
        sum(
            (fs[(op[a], dp[a])] + p[a]) * y_stop[(a, m, n)]
            for (a, m, n) in stop_keys
        )
        -
        sum(p[a] for a in A)
    )

    # Each passenger group assigned at most once
    @constraint(model, [a in A],
        sum((y_dir[(aa, m, n)] for (aa, m, n) in direct_keys if aa == a); init = 0)
        +
        sum((y_stop[(aa, m, n)] for (aa, m, n) in stop_keys if aa == a); init = 0)
        <= 1
    )

    # Capacity on each scheduled leg
    for leg in scheduled
        m = leg.leg_index
        n = leg.plane

        @constraint(model,
            sum(
                (q[a] * y_dir[(a, mm, nn)]
                 for (a, mm, nn) in direct_keys
                 if mm == m && nn == n);
                init = 0
            )
            +
            sum(
                (q[a] * y_stop[(a, mm, nn)]
                 for (a, mm, nn) in stop_keys
                 if (mm == m && nn == n) || (mm + 1 == m && nn == n));
                init = 0
            )
            <= cap_u
        )
    end

    # Direct-only passengers ride alone
    for (a, m, n) in direct_keys
        if so[a] == 0
            others_on_same_leg =
                sum(
                    (y_dir[(aa, mm, nn)]
                     for (aa, mm, nn) in direct_keys
                     if !(aa == a && mm == m && nn == n) && mm == m && nn == n);
                    init = 0
                )
                +
                sum(
                    (y_stop[(aa, mm, nn)]
                     for (aa, mm, nn) in stop_keys
                     if (mm == m && nn == n) || (mm + 1 == m && nn == n));
                    init = 0
                )

            bigM = length(direct_keys) + length(stop_keys)
            @constraint(model, others_on_same_leg <= bigM * (1 - y_dir[(a, m, n)]))
        end
    end

    optimize!(model)

    model.ext[:direct_keys] = direct_keys
    model.ext[:stop_keys] = stop_keys

    return model, scheduled
end

function extract_assignments2(model::Model)
    assignments = PassengerAssignment[]

    direct_keys = model.ext[:direct_keys]
    stop_keys = model.ext[:stop_keys]

    for key in direct_keys
        if value(model[:y_dir][key]) > 0.5
            a, m, n = key
            push!(assignments, PassengerAssignment(a, n, [m]))
        end
    end

    for key in stop_keys
        if value(model[:y_stop][key]) > 0.5
            a, m, n = key
            push!(assignments, PassengerAssignment(a, n, [m, m + 1]))
        end
    end

    return assignments
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

function fitnessFunction(
    evtols::allPlaneSolution,
    assignments::Vector{PassengerAssignment},
    bmax::Float32,
    bmin::Float32,
    dist::Dict{Tuple{Int,Int},Float64},
    ec::Float32,
    battery_per_km::Float32,
    rt::Matrix{Int},
    ET::Int,
    T::Int,
    V::Int,
    cap_flt::Int,
    cap_v::Dict{},
    data
)
    P = FeasibilityCheck(
        bmax,
        bmin,
        dist,
        ec,
        battery_per_km,
        evtols,
        rt,
        ET,
        T,
        V,
        cap_flt,
        cap_v
    )

    A  = data.A
    op = data.op
    dp = data.dp
    fd = data.fd
    fs = data.fs
    c  = data.c
    so = data.so
    p  = data.p

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
            fitnessvalue -= 400.0
        end
    end

    # 5. Large infeasibility penalty
    fitnessvalue -= sum(P) * 1_000_000.0

    return fitnessvalue
end

###############################################################################
# Usage
###############################################################################
# excel_file = joinpath(@__DIR__, "inputData.xlsx")
# data = load_data(excel_file)

# Vmax = maximum(data.V)
# rt = zeros(Int, Vmax, Vmax)
# for i in data.V, j in data.V
#     rt[i, j] = data.rt[(i, j)]
# end

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

# start_time = time()

# best_solutions = generate_best_initial_solutions(data, rt; n_runs=10000, top_k=1, maxLegs=6, maxTurnaround=20)

# elapsed_time = time() - start_time

# println("Time used: ", round(elapsed_time, digits=3), " seconds")

# for (rank, sol) in enumerate(best_solutions)
#     println("====================================")
#     println("Rank: ", rank)
#     println("Run: ", sol.run)
#     println("Fitness: ", sol.fitness)
#     println("P: ", sol.P)
#     println()

#     print_chromosome_table(sol.evtols)
#     println()
#     print_assignments(sol.assignments, data)
#     println()
# end

# ---------------------------------------------- #

function Possible_TurnaroundTime(plane::planeSolution, plane_idx::Int, schedule::Vector{ScheduledLeg}, te::Int, maxTurnaround::Int, ET::Int)

    m = plane.flightLegs
    
    tft = 0

    for leg in schedule
        if leg.plane == plane_idx && leg.leg_index == m
            tft = leg.arr
            break
        end
    end

    # println(tft)
    maxtime = minimum([maximum([ET-tft,plane.turnaroundTime[m]]), maxTurnaround])

    return range(te, maxtime)
end


function Change(plane::planeSolution, new_time::Int, turnaround_idx::Int)

    plane_changed = plane

    plane_changed.turnaroundTime[turnaround_idx] = new_time

    return plane_changed

end

function Best_Change(planes::allPlaneSolution, te::Float64, maxTurnaround::Int, ET::Float64, data, rt)

    temp_sol = 0

    assignment, scheduled = assign_passengersV2(planes, data, rt)

    best_ass = assignment
    best_sol = deepcopy(planes)
    best_obj = fitnessFunction(planes, assignment, Float32(data.bmax), Float32(data.bmin), data.dist, Float32(data.ec), Float32(data.battery_per_km), rt, Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V), Int(round(data.cap_flt)), data.cap_v, data)

    org_obj = copy(best_obj)


    for i in 1:length(best_sol.planes)

        m = best_sol.planes[i].flightLegs

        if m > 1

            turnaroundTimes = Possible_TurnaroundTime(best_sol.planes[i], i, scheduled, Int(te), maxTurnaround, Int(ET))
               
            for j in 1:m

                for t in turnaroundTimes
                    temp_sol = deepcopy(planes)
                    temp_sol.planes[i] = Change(temp_sol.planes[i], t, j)

                    assignment, scheduled = assign_passengersV2(temp_sol, data, rt)
                    temp_obj = fitnessFunction(temp_sol, assignment,  Float32(data.bmax), Float32(data.bmin), data.dist, Float32(data.ec), Float32(data.battery_per_km), rt, Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V), Int(round(data.cap_flt)), data.cap_v, data)

                    if temp_obj > best_obj
                        best_ass = copy(assignment)
                        best_sol = temp_sol
                        best_obj = temp_obj
                    end
                end
            end
        end

       
    end

    println("\nOriginal obj value: $(org_obj)")
    println("\nNew obj value: $(best_obj)")


    return best_sol, best_ass, best_obj

end

###############################################################################
# Simple test for Best_Change function
###############################################################################
function test_best_change()
    println("\nTesting Best_Change function...")
    println("-" ^ 60)
    
    # # Create a simple test solution
    # plane1 = planeSolution(Int32(2), Int32[1, 2, 1], Int32[5, 8])
    # plane2 = planeSolution(Int32(2), Int32[1, 3, 1], Int32[10, 8])
    # test_solution = allPlaneSolution([plane1, plane2])
    
    # # println("Initial solution:")
    # # print_chromosome_table(test_solution)
    
    # # Create a simple scheduled legs list for testing
    # test_schedule = ScheduledLeg[]
    # push!(test_schedule, ScheduledLeg(1, 1, 1, 2, 5, 20, 5))
    # push!(test_schedule, ScheduledLeg(1, 2, 2, 1, 28, 43, 5))

    # push!(test_schedule, ScheduledLeg(2, 1, 1, 3, 10, 17, 5))
    # push!(test_schedule, ScheduledLeg(2, 2, 3, 1, 25, 32, 5))

    
    # te = 5
    # maxTurnaround = 20
    # ET = 120

    excel_file = joinpath(@__DIR__, "inputData.xlsx")
    data = load_data(excel_file)

    Vmax = maximum(data.V)
    rt = zeros(Int, Vmax, Vmax)
    for i in data.V, j in data.V
        rt[i, j] = data.rt[(i, j)]
    end


    maxTurnaround=30

    start_time = time()
    Best_sols = generate_best_initial_solutions(data, rt, top_k = 1)
    elapsed_time = time() - start_time
    println("Run 1000 inits solutions time: $(round(elapsed_time, digits=4)) seconds")



    evtols_init = Best_sols[1].evtols
    assignment_init = Best_sols[1].assignments

    println("\nTesting Best Change function")
    start_time = time()
    sol_out, ass_out, obj_out = Best_Change(evtols_init, data.te, maxTurnaround, data.ET, data, rt)
    elapsed_time = time() - start_time
    println("Best_Change execution time: $(round(elapsed_time, digits=4)) seconds")


    println("Initial solution:")
    print_chromosome_table(evtols_init)
    print_assignments(assignment_init, data)


    println("New solution:")
    print_chromosome_table(sol_out)
    print_assignments(ass_out, data)


    println("-" ^ 60)
end

# Uncomment to run the test:
# test_best_change()



function Destructor(plane::planeSolution, idxs::Vector{Int64})
    m = Int(plane.flightLegs)

    if m == 2 && !isempty(idxs)
        base = plane.route[1]
        resize!(plane.route, 1)
        plane.route[1] = base
        empty!(plane.turnaroundTime)
        plane.flightLegs = Int32(0)
        return
    end

    del = sort(unique(Int.(idxs)))

    # keep only internal stop indices for safety
    del = [idx for idx in del if 2 <= idx <= m]
    isempty(del) && return

    plane.flightLegs -= Int32(length(del))
    deleteat!(plane.route, del)
    deleteat!(plane.turnaroundTime, del)
end

function has_consecutive_duplicates(route)
    for i in 1:(length(route)-1)
        if route[i] == route[i+1]
            return true
        end
    end
    return false
end



function Constructor(plane::planeSolution, VP::Int, idx::Int, data)
    te32 = Int32(data.te)

    if plane.flightLegs == 0
        base = Int32(plane.route[1])

        # Build first tour as base -> VP -> base
        insert!(plane.route, 2, Int32(VP))
        insert!(plane.route, 3, base)

        insert!(plane.turnaroundTime, 1, te32)
        insert!(plane.turnaroundTime, 2, te32)

        plane.flightLegs += 2
        return
    end

    plane.flightLegs += 1
    insert!(plane.route, idx, Int32(VP))
    insert!(plane.turnaroundTime, idx, te32)

end

function two_opt(plane::planeSolution, idx1::Int, idx2::Int)

    # Nothing to improve if there are fewer than 2 flown legs.
    if plane.flightLegs < 2
        return plane
    end

    # Ensure ordered bounds.
    if idx1 > idx2
        idx1, idx2 = idx2, idx1
    end

    # Keep first/last node fixed (base), reverse only internal part.
    lo = max(2, idx1)
    hi = min(length(plane.route) - 1, idx2)

    if lo >= hi
        return plane
    end

    rev_seg = reverse(plane.route[lo:hi])

    # Reject move if it creates consecutive identical vertiports.
    if lo > 1 && plane.route[lo - 1] == rev_seg[1]
        return plane
    end
    if hi < length(plane.route) && rev_seg[end] == plane.route[hi + 1]
        return plane
    end
    for k in 1:(length(rev_seg) - 1)
        if rev_seg[k] == rev_seg[k + 1]
            return plane
        end
    end

    plane.route[lo:hi] = rev_seg
    return plane
end

function UpdateTurnAroundTimes(planes::allPlaneSolution, from::Int64, maxTurnaround::Int64, data)
    op = data.op
    te = Int(data.te)
    dt = data.dt
    rt = data.rt

    for plane in planes.planes
        flightLegs = Int(plane.flightLegs)
        if flightLegs == 0
            continue
        end

        route = plane.route
        from_k = clamp(Int(from), 1, flightLegs)

        # Rebuild state up to the first leg we want to update
        current_time = 0
        current_VP = Int(route[1])

        for k in 1:(from_k - 1)
            next_VP = Int(route[k + 1])
            current_time += Int(plane.turnaroundTime[k]) + rt[(current_VP, next_VP)]
            current_VP = next_VP
        end

        # Recompute turnaround times from from_k onward
        for k in from_k:flightLegs
            candidate_pass = [
                dt[a] for (a, v) in op
                if v == current_VP && current_time + te <= dt[a] <= current_time + maxTurnaround
            ]

            if !isempty(candidate_pass)
                # deterministic: latest feasible passenger time in the window
                chosen_time = maximum(candidate_pass) - current_time
                plane.turnaroundTime[k] = Int32(round(Int, chosen_time))
            else
                plane.turnaroundTime[k] = Int32(rand(te:maxTurnaround))
            end

            next_VP = Int(route[k + 1])
            current_time += Int(plane.turnaroundTime[k]) + rt[(current_VP, next_VP)]
            current_VP = next_VP
        end
    end
end

function obj(planes::allPlaneSolution, data, rt)
    assignments, scheduled = assign_passengersV2(planes, data, rt)

    fitness = fitnessFunction(
            planes,
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

    return fitness
end

function DestructLoop(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    N = data.N
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Do not delete first node (route[1]); last route node is not in this range anyway
        for idx in 2:m
            temp_sol = deepcopy(planes)

            Destructor(temp_sol.planes[n], [idx])

            if has_consecutive_duplicates(temp_sol.planes[n].route)
                continue
            end

            new_obj = obj(temp_sol, data, rt)

            temp_sol2 = deepcopy(temp_sol)
            from_idx = Int64(max(1, idx - 1))
            UpdateTurnAroundTimes(temp_sol2, from_idx, maxTurnaround, data)
            new_obj2 = obj(temp_sol2, data, rt)

            if new_obj2 > new_obj
                temp_sol = temp_sol2
                new_obj = new_obj2
            end

            if new_obj > best_obj
                best_obj = new_obj
                best_sol = deepcopy(temp_sol)
            end
        end
    end

    return best_obj, best_sol
end

function ConstructLoop(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)

    N = data.N
    V = data.V
    best_obj = init_obj
    best_sol = deepcopy(planes)


    for n in N
        m = Int(planes.planes[n].flightLegs)

        if m == 0
            base = Int(planes.planes[n].route[1])

            for v in V
                if v == base
                    continue
                end

                temp_sol = deepcopy(planes)
                Constructor(temp_sol.planes[n], v, 2, data)

                new_obj = obj(temp_sol, data, rt)

                temp_sol2 = deepcopy(temp_sol)
                UpdateTurnAroundTimes(temp_sol2, 1, maxTurnaround, data)
                new_obj2 = obj(temp_sol2, data, rt)

                if new_obj2 > new_obj
                    temp_sol = temp_sol2
                    new_obj = new_obj2
                end

                if new_obj > best_obj
                    best_obj = new_obj
                    best_sol = deepcopy(temp_sol)
                end
            end

            continue
        end

        for idx in 2:(m + 1)
            for v in V
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

                    if new_obj > best_obj
                        best_obj = new_obj
                        best_sol = deepcopy(temp_sol)
                    end
                end
            end
        end
    end

    return best_obj, best_sol
end

function Swap(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)
    N = data.N
    V = data.V
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Need at least 2 legs to have an internal stop to swap
        if m <= 1
            continue
        end

        # Internal route positions only (avoid route[1] base)
        for idx in 2:m
            # Neighbors around the original stop at idx
            left_v = planes.planes[n].route[idx - 1]
            old_v = planes.planes[n].route[idx]
            right_v = planes.planes[n].route[idx + 1]

            for v in V
                # Skip no-op and adjacent duplicates
                if v == old_v || v == left_v || v == right_v
                    continue
                end

                cand = deepcopy(planes)

                # 1) Insert new stop at idx
                Constructor(cand.planes[n], v, idx, data)

                # 2) Remove original stop (shifted to idx+1 after insertion)
                Destructor(cand.planes[n], [Int64(idx + 1)])

                new_obj = obj(cand, data, rt)

                # Candidate with updated turnaround times
                cand2 = deepcopy(cand)
                from_idx = Int64(max(1, idx - 1))
                UpdateTurnAroundTimes(cand2, from_idx, maxTurnaround, data)
                new_obj2 = obj(cand2, data, rt)

                if new_obj2 > new_obj
                    cand = cand2
                    new_obj = new_obj2
                end

                if new_obj > best_obj
                    best_obj = new_obj
                    best_sol = deepcopy(cand)
                end
            end
        end
    end

    return best_obj, best_sol
end

function two_opt_Loop(planes::allPlaneSolution, maxTurnaround::Int64, init_obj::Float64, data, rt)
    N = data.N
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Need at least 3 legs to have two internal indices to reverse.
        if m <= 2
            continue
        end

        for i in 2:(m - 1)
            for j in (i + 1):m
                temp_sol = deepcopy(planes)
                two_opt(temp_sol.planes[n], i, j)

                temp_obj = obj(temp_sol, data, rt)

                # Re-evaluate with refreshed turnaround times from the changed area.
                temp_sol2 = deepcopy(temp_sol)
                from_idx = Int64(max(1, i - 1))
                UpdateTurnAroundTimes(temp_sol2, from_idx, maxTurnaround, data)
                temp_obj2 = obj(temp_sol2, data, rt)

                if temp_obj2 > temp_obj
                    temp_sol = temp_sol2
                    temp_obj = temp_obj2
                end

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                end
            end
        end
    end
    
    return best_obj, best_sol
end


function Heuristic(maxTurnaround::Int64, MaxTime::Int32, data, rt)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    while elapsed <= Float64(MaxTime)
        nr = 1
        Best_sols = generate_best_initial_solutions(data, rt, top_k = nr, n_runs = 100, maxTurnaround = maxTurnaround)

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

excel_file = joinpath(@__DIR__, "inputData.xlsx")
data = load_data(excel_file)

Vmax = maximum(data.V)
rt = zeros(Int, Vmax, Vmax)
for i in data.V, j in data.V
    rt[i, j] = data.rt[(i, j)]
end


maxTurnaround = 50
Maxtime = Int32(10) 
nr = 1

(best_obj, best_sol, iterations) = Heuristic(maxTurnaround, Maxtime, data, rt)

println("Heuristic ran $(iterations) iterations")
println("Best solution:")
println("Objective Value: $(best_obj)")
print_chromosome_table(best_sol)

assignments, scheduled = assign_passengersV2(best_sol, data, rt)

# println("Passenger Assignment")
print_assignments(assignments, data)

print_schedule_pretty(scheduled)
# println("Max waiting time: $(data.w)")

# start_time = time()
# Best_sols = generate_best_initial_solutions(data, rt, top_k = nr)

# Org_obj = Best_sols[nr].fitness

# evtols_init = Best_sols[nr].evtols


# println("Initial solution:"),
# print_chromosome_table(evtols_init)

# # (best_obj, best_sol) = DestructLoop(evtols_init, maxTurnaround, Org_obj, data, rt)
# (best_obj, best_sol) = ConstructLoop(evtols_init, maxTurnaround, Org_obj, data, rt)
# # (best_obj, best_sol) = Swap(evtols_init, maxTurnaround, Org_obj, data, rt)

# println("Best solution:")
# print_chromosome_table(best_sol)

# println("Initial Obj $(Org_obj)")
# println("Best Obj $(best_obj)")