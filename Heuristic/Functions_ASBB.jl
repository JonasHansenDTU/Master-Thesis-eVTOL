using XLSX
using DataFrames
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
    infra = read_sheet(excel_file, "Infrastructure")
    pax   = read_sheet(excel_file, "PassengerGroups")
    plane = read_sheet_any(excel_file, ["PlaneData"])

    ###########################################################################
    # Infrastructure columns
    ###########################################################################
    id_col     = find_col(infra, [:id])
    pads_col   = find_col(infra, [:number_of_parking_pads, :parking_pads, :pads])

    # Coordinates can be either one string column "coordinates"
    # or two numeric columns such as "latitude", "longitude".
    coord_col = if any(Symbol(String(n)) == :coordinates for n in names(infra))
        find_col(infra, [:coordinates])
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

    weights = build_vertiport_weights(V, op, dp, A)

    planes = planeSolution[]

    for n in N
        base = bv[n]

        # choose number of legs
        allowed_legs = [i for i in 0:maxLegs if i != 1]
        flightLegs = rand(allowed_legs)

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
        for k in 1:flightLegs
            push!(turnaroundTime, Int32(rand(te:maxTurnaround)))
        end

        push!(planes, planeSolution(
            Int32(flightLegs),
            route,
            turnaroundTime
        ))
    end

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

        for leg2 in scheduled
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
        if a in assigned_groups
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
excel_file = joinpath(@__DIR__, "inputData.xlsx")
data = load_data(excel_file)

Vmax = maximum(data.V)
rt = zeros(Int, Vmax, Vmax)
for i in data.V, j in data.V
    rt[i, j] = data.rt[(i, j)]
end

function generate_best_initial_solutions(data, rt; n_runs::Int=1000, top_k::Int=10, maxLegs::Int=5, maxTurnaround::Int=30, print_each::Bool=false)
    results = NamedTuple[]

    for run in 1:n_runs
        evtols_init = initial_chromosome_solution(data; maxLegs=maxLegs, maxTurnaround=maxTurnaround)

        assignments, scheduled = assign_passengers(evtols_init, data, rt)

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

start_time = time()

best_solutions = generate_best_initial_solutions(data, rt; n_runs=10000, top_k=1, maxLegs=6, maxTurnaround=20)

elapsed_time = time() - start_time

println("Time used: ", round(elapsed_time, digits=3), " seconds")

for (rank, sol) in enumerate(best_solutions)
    println("====================================")
    println("Rank: ", rank)
    println("Run: ", sol.run)
    println("Fitness: ", sol.fitness)
    println("P: ", sol.P)
    println()

    print_chromosome_table(sol.evtols)
    println()
    print_assignments(sol.assignments, data)
    println()
end


function Change(planes::allPlaneSolution, te::Int, maxTurnaround::Int, schedule::Vector{ScheduledLeg}, ET::Int)

    plane_idx = rand(1:length(planes.planes))
    plane = planes.planes[plane_idx]

    turnaround_idx = rand(1:length(plane.turnaroundTime))

    m = plane.flightLegs
    
    tft = 0

    for leg in schedule
        if leg.plane == plane_idx && leg.leg_index == m
            tft = leg.arr
            break
        end
    end

    println(tft)
    maxtime = maximum([ET-tft,plane.turnaroundTime[m]])

    t = rand(te:maxtime)

    plane.turnaroundTime[turnaround_idx] = t

    println(te)
    println(t)
    println(maxtime)
end


function insert(plane::allPlaneSolution, data; maxTurnaround=30)
    planeidx = rand(1:length(plane.planes))
    p = plane.planes[planeidx]

    route = p.route

    # Need at least 2 nodes to insert between
    if length(route) < 2
        return false
    end

    # Choose insertion index (not first or last)
    idx = rand(2:length(route))

    prev = Int(route[idx - 1])
    next = Int(route[idx])

    # Exclude neighbors
    forbidden = [prev, next]
    candidates = [v for v in data.V if !(v in forbidden)]

    if isempty(candidates)
        return false
    end

    newvertiport = rand(candidates)
    newturnaround = rand(Int(round(data.te)):maxTurnaround)

    insert!(p.route, idx, Int32(newvertiport))
    insert!(p.turnaroundTime, idx - 1, Int32(newturnaround))

    p.flightLegs += 1

    return true
end

println("\n===== TEST INSERT FUNCTION =====")

# Create one random solution
evtols = deepcopy(best_solutions[1].evtols)
println("\n--- BEFORE ---")
print_chromosome_table(evtols)

# Copy for comparison (important!)
evtols_before = deepcopy(evtols)

# Apply insert
success = insert(evtols, data)

println("\nInsert success: ", success)

println("\n--- AFTER ---")
print_chromosome_table(evtols)

function delete(population::allPlaneSolution, data)

    # Only planes with removable intermediate nodes, route length > 2
    candidates = [
        i for i in 1:length(population.planes)
        if length(population.planes[i].route) > 2
    ]

    if isempty(candidates)
        return false
    end

    # Try until we find a valid deletion
    for _ in 1:10  # small retry loop to avoid infinite failure
        pidx = rand(candidates)
        p = population.planes[pidx]

        # do not dele in the first or last leg 
        idx = rand(2:length(p.route)-1)

        prev = p.route[idx - 1]
        next = p.route[idx + 1]

        # Do not delete if it creates a degenerate edge (prev -> same -> next)
        if prev == next
            continue
        end

        deleteat!(p.route, idx)
        deleteat!(p.turnaroundTime, idx - 1)

        # updatte flight legs count
        p.flightLegs = length(p.route) - 1

        return true
    end

    return false
end

println("\n===== TEST DELETE FUNCTION =====")

evtols = deepcopy(best_solutions[1].evtols)

println("\n--- BEFORE ---")
print_chromosome_table(evtols)

delete(evtols, data)

println("\n--- AFTER ---")
print_chromosome_table(evtols)

for p in evtols.planes
    println("route length = ", length(p.route),
            " | legs = ", p.flightLegs,
            " | turnaround = ", length(p.turnaroundTime))
end
