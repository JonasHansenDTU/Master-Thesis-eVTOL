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
        cap_v = cap_v, cap_u = cap_u,
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

function FeasibilityCheck(bmax::Float32, bmin::Float32,
    dist::Dict{Tuple{Int,Int},Float64},
    ec::Float32, battery_per_km::Float32,
    evtols::allPlaneSolution,
    rt::Matrix{Int}, ET::Int, T::Int, V::Int, cap_v::Dict{})

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

    for n in 1:nPlanes-1
        allowed_legs[n+1] = (rand(Set_of_flightlegs))
        
        while n > 1 && (allowed_legs[n+1] == allowed_legs[n] -1 
                    || allowed_legs[n+1] == allowed_legs[n] +1) 
            allowed_legs[n+1] = rand(Set_of_flightlegs) 
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
                chosen_time = Candidate_Pass[end] - current_time
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
    if so == 1
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
    Price_sort = sort(A, by = a -> price_by_group[a], rev = true)  # descending

    for a in Price_sort
        if so[a] == 0
            ass = find_direct_leg!V2(scheduled, a, op, dp, dt, q, w, so)
            if ass !== nothing
                push!(assignments, ass)
                push!(assigned_groups, a)
            end
        end
        if so[a] == 1
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

function fitnessFunction(evtols::allPlaneSolution,assignments::Vector{PassengerAssignment},bmax::Float32,bmin::Float32,dist::Dict{Tuple{Int,Int},Float64},ec::Float32,
                        battery_per_km::Float32, rt::Matrix{Int}, ET::Int, T::Int, V::Int, cap_v::Dict{}, data)
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

function initial_chromosome_solution2(data, rt; maxLegs::Int=5, maxTurnaround::Int=30, top_c::Int=3)
    V  = data.V
    N  = data.N
    A  = data.A
    bv = data.bv
    op = data.op
    dp = data.dp
    dt = data.dt
    q  = data.q
    so = data.so
    w  = Int(round(data.w))
    te = Int(round(data.te))
    fd = data.fd
    fs = data.fs
    cap_u = Int(round(data.cap_u))

    weights = build_vertiport_weights(V, op, dp, A)

    planes = Vector{Union{Nothing,planeSolution}}(undef, length(N))
    assignments = PassengerAssignment[]
    assigned_groups = Set{Int}()

    # ----------------------------
    # Random construction mode
    # ----------------------------
    early_first = rand() < 0.5

    passenger_value = Dict(a => (so[a] == 1 ? fs[(op[a], dp[a])] : fd[(op[a], dp[a])]) for a in A)

    if early_first
        sorted_groups = sort(A, by = a -> (dt[a], -passenger_value[a]))
    else
        sorted_groups = sort(A, by = a -> (-dt[a], -passenger_value[a]))
    end

    # ----------------------------
    # Random plane order
    # ----------------------------
    plane_order = shuffle(collect(1:length(N)))
    nPlanes = length(N)

    # ----------------------------
    # Distribute total flight legs
    # ----------------------------
    choices = [i for i in 0:(maxLegs * nPlanes) if i != 1]
    total_flight_legs = rand(choices)

    allowed_legs = zeros(Int, nPlanes)

    while true
        remaining = total_flight_legs
        allowed_legs .= 0
        feasible_distribution = true

        for n in 1:(nPlanes - 1)
            feasible_choices = Int[]

            for x in 0:maxLegs
                rem_after = remaining - x
                if x != 1 && rem_after >= 0
                    max_possible_rest = maxLegs * (nPlanes - n)
                    if (rem_after == 0 || rem_after >= 2) && rem_after <= max_possible_rest
                        push!(feasible_choices, x)
                    end
                end
            end

            if isempty(feasible_choices)
                feasible_distribution = false
                break
            end

            choice = rand(feasible_choices)
            allowed_legs[n] = choice
            remaining -= choice
        end

        if feasible_distribution &&
           remaining <= maxLegs &&
           remaining != 1 &&
           remaining >= 0
            allowed_legs[nPlanes] = remaining
            allowed_legs = shuffle!(allowed_legs)
            break
        end
    end

    # ----------------------------
    # Phase 1: build routes + seed assignments
    # ----------------------------
    for n in plane_order
        base = bv[n]
        flightLegs_target = allowed_legs[n]

        route = Int32[base]
        turnaroundTime = Int32[]

        current_vp = base
        current_time = 0
        legs_used = 0

        while legs_used < flightLegs_target
            candidates = NamedTuple[]

            for a in sorted_groups
                if a in assigned_groups
                    continue
                end

                # ------------------------
                # Direct option
                # ------------------------
                if current_vp == op[a]
                    if legs_used + 1 <= flightLegs_target
                        wait1 = max(te, Int(round(dt[a] - current_time)))
                        wait1 = min(wait1, maxTurnaround)
                        dep1 = current_time + wait1

                        if dt[a] <= dep1 <= dt[a] + w
                            score = passenger_value[a] - 0.1 * wait1
                            push!(candidates, (
                                score = score,
                                group = a,
                                mode = :direct,
                                midpoint = nothing,
                                reposition = false,
                                rep_wait = 0,
                                wait1 = wait1,
                                wait2 = 0
                            ))
                        end
                    end
                else
                    if legs_used + 2 <= flightLegs_target
                        rep_wait = te
                        arrival_origin = current_time + rep_wait + rt[current_vp, op[a]]
                        wait1 = max(te, Int(round(dt[a] - arrival_origin)))
                        wait1 = min(wait1, maxTurnaround)
                        dep1 = arrival_origin + wait1

                        if dt[a] <= dep1 <= dt[a] + w
                            score = passenger_value[a] - 0.1 * wait1 - 0.2 * rt[current_vp, op[a]]
                            push!(candidates, (
                                score = score,
                                group = a,
                                mode = :direct,
                                midpoint = nothing,
                                reposition = true,
                                rep_wait = rep_wait,
                                wait1 = wait1,
                                wait2 = 0
                            ))
                        end
                    end
                end

                # ------------------------
                # Stopover option
                # ------------------------
                if so[a] == 1
                    mids = [v for v in V if v != op[a] && v != dp[a]]
                    if !isempty(mids)
                        chosen_mid = weighted_choice(mids, [weights[findfirst(==(v), V)] for v in mids])

                        if current_vp == op[a]
                            if legs_used + 2 <= flightLegs_target
                                wait1 = max(te, Int(round(dt[a] - current_time)))
                                wait1 = min(wait1, maxTurnaround)
                                dep1 = current_time + wait1
                                arr1 = dep1 + rt[op[a], chosen_mid]
                                wait2 = te
                                dep2 = arr1 + wait2

                                if dt[a] <= dep1 <= dt[a] + w && dep2 >= arr1
                                    score = passenger_value[a] - 0.1 * (wait1 + wait2)
                                    push!(candidates, (
                                        score = score,
                                        group = a,
                                        mode = :stop,
                                        midpoint = chosen_mid,
                                        reposition = false,
                                        rep_wait = 0,
                                        wait1 = wait1,
                                        wait2 = wait2
                                    ))
                                end
                            end
                        else
                            if legs_used + 3 <= flightLegs_target
                                rep_wait = te
                                arrival_origin = current_time + rep_wait + rt[current_vp, op[a]]
                                wait1 = max(te, Int(round(dt[a] - arrival_origin)))
                                wait1 = min(wait1, maxTurnaround)
                                dep1 = arrival_origin + wait1
                                arr1 = dep1 + rt[op[a], chosen_mid]
                                wait2 = te
                                dep2 = arr1 + wait2

                                if dt[a] <= dep1 <= dt[a] + w && dep2 >= arr1
                                    score = passenger_value[a] - 0.1 * (wait1 + wait2) - 0.2 * rt[current_vp, op[a]]
                                    push!(candidates, (
                                        score = score,
                                        group = a,
                                        mode = :stop,
                                        midpoint = chosen_mid,
                                        reposition = true,
                                        rep_wait = rep_wait,
                                        wait1 = wait1,
                                        wait2 = wait2
                                    ))
                                end
                            end
                        end
                    end
                end
            end

            if isempty(candidates)
                break
            end

            candidates = sort(candidates, by = x -> x.score, rev = true)
            pool = candidates[1:min(top_c, length(candidates))]
            chosen = rand(pool)

            a = chosen.group

            if chosen.reposition
                push!(route, Int32(op[a]))
                push!(turnaroundTime, Int32(chosen.rep_wait))
                current_time += turnaroundTime[end] + rt[current_vp, op[a]]
                current_vp = op[a]
                legs_used += 1
            end

            if chosen.mode == :direct
                service_leg = legs_used + 1

                push!(route, Int32(dp[a]))
                push!(turnaroundTime, Int32(chosen.wait1))
                current_time += turnaroundTime[end] + rt[current_vp, dp[a]]
                current_vp = dp[a]
                legs_used += 1

                push!(assignments, PassengerAssignment(a, n, [service_leg]))
                push!(assigned_groups, a)

            elseif chosen.mode == :stop
                leg1 = legs_used + 1
                leg2 = legs_used + 2

                push!(route, Int32(chosen.midpoint))
                push!(turnaroundTime, Int32(chosen.wait1))
                current_time += turnaroundTime[end] + rt[current_vp, chosen.midpoint]
                current_vp = chosen.midpoint
                legs_used += 1

                push!(route, Int32(dp[a]))
                push!(turnaroundTime, Int32(chosen.wait2))
                current_time += turnaroundTime[end] + rt[current_vp, dp[a]]
                current_vp = dp[a]
                legs_used += 1

                push!(assignments, PassengerAssignment(a, n, [leg1, leg2]))
                push!(assigned_groups, a)
            end
        end

        # If no more passenger candidate: either fill spare legs randomly or return to base
        while legs_used < flightLegs_target - 1
            candidates = [v for v in V if v != current_vp]
            if isempty(candidates)
                break
            end

            nxt = weighted_choice(candidates, [weights[findfirst(==(v), V)] for v in candidates])
            push!(route, Int32(nxt))
            push!(turnaroundTime, Int32(rand(te:maxTurnaround)))
            current_time += turnaroundTime[end] + rt[current_vp, nxt]
            current_vp = nxt
            legs_used += 1
        end

        if legs_used < flightLegs_target && current_vp != base
            push!(route, Int32(base))
            push!(turnaroundTime, Int32(rand(te:maxTurnaround)))
            current_time += turnaroundTime[end] + rt[current_vp, base]
            current_vp = base
            legs_used += 1
        end

        planes[n] = planeSolution(
            Int32(length(route) - 1),
            route,
            turnaroundTime
        )
    end

    evtols_init = allPlaneSolution(planeSolution[p for p in planes])
    scheduled = build_scheduled_legs(evtols_init, rt, cap_u)

    # ----------------------------
    # Phase 2: assign remaining groups on built schedule
    # ----------------------------
    assigned_groups = Set(ass.group for ass in assignments)

    remaining_capacity = Dict((leg.plane, leg.leg_index) => cap_u for leg in scheduled)
    exclusive_leg = Dict((leg.plane, leg.leg_index) => false for leg in scheduled)

    for ass in assignments
        a = ass.group
        for legidx in ass.legs
            remaining_capacity[(ass.plane, legidx)] -= q[a]
        end
        if so[a] == 0
            @assert length(ass.legs) == 1
            exclusive_leg[(ass.plane, ass.legs[1])] = true
        end
    end

    for a in sorted_groups
        if a in assigned_groups
            continue
        end

        earliest = Int(round(dt[a]))
        latest = earliest + w

        if so[a] == 0
            chosen = nothing
            best_arr = typemax(Int)

            for leg in scheduled
                key = (leg.plane, leg.leg_index)
                if leg.from == op[a] &&
                   leg.to == dp[a] &&
                   earliest <= leg.dep <= latest &&
                   remaining_capacity[key] == cap_u &&
                   !exclusive_leg[key]

                    if leg.arr < best_arr
                        best_arr = leg.arr
                        chosen = leg
                    end
                end
            end

            if chosen !== nothing
                key = (chosen.plane, chosen.leg_index)
                remaining_capacity[key] -= q[a]
                exclusive_leg[key] = true
                push!(assignments, PassengerAssignment(a, chosen.plane, [chosen.leg_index]))
                push!(assigned_groups, a)
            end
        else
            chosen_direct = nothing
            best_direct_arr = typemax(Int)

            for leg in scheduled
                key = (leg.plane, leg.leg_index)
                if leg.from == op[a] &&
                   leg.to == dp[a] &&
                   earliest <= leg.dep <= latest &&
                   remaining_capacity[key] >= q[a] &&
                   !exclusive_leg[key]

                    if leg.arr < best_direct_arr
                        best_direct_arr = leg.arr
                        chosen_direct = leg
                    end
                end
            end

            if chosen_direct !== nothing
                key = (chosen_direct.plane, chosen_direct.leg_index)
                remaining_capacity[key] -= q[a]
                push!(assignments, PassengerAssignment(a, chosen_direct.plane, [chosen_direct.leg_index]))
                push!(assigned_groups, a)
                continue
            end

            chosen_pair = nothing
            best_pair_arr = typemax(Int)

            for leg1 in scheduled
                key1 = (leg1.plane, leg1.leg_index)
                if leg1.from != op[a]
                    continue
                end
                if !(earliest <= leg1.dep <= latest)
                    continue
                end
                if remaining_capacity[key1] < q[a] || exclusive_leg[key1]
                    continue
                end

                for leg2 in scheduled
                    key2 = (leg2.plane, leg2.leg_index)
                    if leg2.plane != leg1.plane
                        continue
                    end
                    if leg2.leg_index != leg1.leg_index + 1
                        continue
                    end
                    if leg2.from != leg1.to
                        continue
                    end
                    if leg2.to != dp[a]
                        continue
                    end
                    if leg2.dep < leg1.arr
                        continue
                    end
                    if remaining_capacity[key2] < q[a] || exclusive_leg[key2]
                        continue
                    end

                    if leg2.arr < best_pair_arr
                        best_pair_arr = leg2.arr
                        chosen_pair = (leg1, leg2)
                    end
                end
            end

            if chosen_pair !== nothing
                leg1, leg2 = chosen_pair
                remaining_capacity[(leg1.plane, leg1.leg_index)] -= q[a]
                remaining_capacity[(leg2.plane, leg2.leg_index)] -= q[a]
                push!(assignments, PassengerAssignment(a, leg1.plane, [leg1.leg_index, leg2.leg_index]))
                push!(assigned_groups, a)
            end
        end
    end

    return evtols_init, assignments, scheduled
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

function generate_best_initial_solutions(data, rt; n_runs::Int=1000, top_k::Int=10, top_c::Int=3,  maxLegs::Int=5, maxTurnaround::Int=30, print_each::Bool=false)
    results = NamedTuple[]

    for run in 1:n_runs
        evtols_init, assignments, scheduled = initial_chromosome_solution2(data, rt; maxLegs=maxLegs, maxTurnaround=maxTurnaround, top_c=top_c)

        #assignments, scheduled = assign_passengersV2(evtols_init, data, rt)

        P = FeasibilityCheck(Float32(data.bmax),Float32(data.bmin),data.dist,Float32(data.ec),Float32(data.battery_per_km),
                    evtols_init,rt,Int(round(data.ET)),maximum(Int.(data.T)),maximum(data.V),data.cap_v)

        fitness = fitnessFunction(evtols_init,assignments,Float32(data.bmax),Float32(data.bmin),data.dist, Float32(data.ec),
                        Float32(data.battery_per_km), rt, Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V), data.cap_v, data)

        push!(results, (run = run, fitness = fitness, evtols = evtols_init, assignments = assignments, scheduled = scheduled, P = P))

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

function Change(plane::planeSolution, new_time::Int, turnaround_idx::Int)

    plane_changed = plane

    plane_changed.turnaroundTime[turnaround_idx] = new_time

    return plane_changed

end

function Destructor(plane::planeSolution, idxs::Vector{Int64}) # Deletes stops (idxs) of one plane
    m = Int(plane.flightLegs)

    # Special case: if route has exactly 2 legs (base -> x -> base),
    # any deletion request should remove the whole tour.
    if m == 2 && !isempty(idxs)
        base = plane.route[1]
        resize!(plane.route, 1)
        plane.route[1] = base
        empty!(plane.turnaroundTime)
        plane.flightLegs = Int32(0)
        return
    end

    # Normal behavior
    del = sort(unique(Int.(idxs)))
    plane.flightLegs -= Int32(length(del))
    deleteat!(plane.route, del)
    deleteat!(plane.turnaroundTime, del)
end

function Constructor(plane::planeSolution, maxTurnaround, VP::Int, idx::Int, data)
    te32 = rand(Int32(data.te):maxTurnaround)

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

    fitness = fitnessFunction(planes,assignments,Float32(data.bmax),Float32(data.bmin),data.dist,Float32(data.ec),
            Float32(data.battery_per_km),rt, Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V), data.cap_v,data)

    return fitness
end

function DestructLoop(planes::allPlaneSolution, init_obj::Float64, data, rt)

    N = data.N
    best_obj = init_obj
    best_sol = deepcopy(planes)

    for n in N
        m = Int(planes.planes[n].flightLegs)

        # Do not delete first node (route[1]); last route node is not in this range anyway
        for idx in 2:m
            temp_sol = deepcopy(planes)

            Destructor(temp_sol.planes[n], [idx])
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

function ConstructLoop(planes::allPlaneSolution, maxTurnaround, init_obj::Float64, data, rt)

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
                Constructor(temp_sol.planes[n],maxTurnaround, v, 2, data)

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
                    Constructor(temp_sol.planes[n],maxTurnaround, v, idx, data)
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
                Constructor(cand.planes[n], maxTurnaround, v, idx, data)

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

function Heuristic(data, rt, maxTurnaround, runtime, top_c)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    while elapsed <= Float64(runtime)
        nr = rand(1:5)
        Best_sols = generate_best_initial_solutions(data, rt; n_runs = 100, top_k = nr, top_c, maxLegs=8, maxTurnaround)

        temp_obj = Best_sols[nr].fitness
        temp_sol = deepcopy(Best_sols[nr].evtols)

        if temp_obj > best_obj
            best_obj = temp_obj
            best_sol = deepcopy(temp_sol)
            println("New best Obj $(best_obj)")
        end

        improvement = true
        while improvement && elapsed <= Float64(runtime)
            improvement = false

            des_obj, des_sol = DestructLoop(temp_sol, temp_obj, data, rt)
            con_obj, con_sol = ConstructLoop(temp_sol, maxTurnaround, temp_obj, data, rt)
            swap_obj, swap_sol = Swap(temp_sol, maxTurnaround, temp_obj, data, rt)

            cand_obj = des_obj
            cand_sol = des_sol

            if swap_obj > cand_obj
                cand_obj = swap_obj
                cand_sol = swap_sol
            end

            if con_obj > cand_obj
                cand_obj = con_obj
                cand_sol = con_sol
            end


            if cand_obj > temp_obj
                temp_obj = cand_obj
                temp_sol = deepcopy(cand_sol)
                improvement = true

                if temp_obj > best_obj
                    best_obj = temp_obj
                    best_sol = deepcopy(temp_sol)
                    println("New best Obj $(best_obj)")
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

function Heuristic2(data, rt, maxTurnaround, runtime, top_c)

    start_ns = time_ns()
    elapsed = 0.0
    iterations = 0

    best_obj = -Inf
    best_sol = allPlaneSolution(planeSolution[])

    while elapsed <= Float64(runtime)
        nr = rand(1:5)
        Best_sols = generate_best_initial_solutions(data, rt; n_runs = 100, top_k = nr, top_c, maxLegs=8, maxTurnaround)

        temp_obj = Best_sols[nr].fitness
        temp_sol = deepcopy(Best_sols[nr].evtols)

        if temp_obj > best_obj
            best_obj = temp_obj
            best_sol = deepcopy(temp_sol)
            println("New best Obj $(best_obj)")
        end

        improvement = true
        while improvement
            elapsed = (time_ns() - start_ns) / 1e9
            if elapsed > Float64(runtime)
                break
            end

            improvement = false

            des_obj, des_sol = DestructLoop(temp_sol, temp_obj, data, rt)
            if des_obj > temp_obj
                temp_obj = des_obj
                temp_sol = des_sol
                improvement = true
            else
                con_obj, con_sol = ConstructLoop(temp_sol, maxTurnaround, temp_obj, data, rt)
                if con_obj > temp_obj
                    temp_obj = con_obj
                    temp_sol = con_sol
                    improvement = true
                else
                    swap_obj, swap_sol = Swap(temp_sol, maxTurnaround, temp_obj, data, rt)
                    if swap_obj > temp_obj
                        temp_obj = swap_obj
                        temp_sol = swap_sol
                        improvement = true
                    end
                end
            end

            if improvement && temp_obj > best_obj
                best_obj = temp_obj
                best_sol = deepcopy(temp_sol)
                println("New best Obj $(best_obj)")
            end
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

nr = 1
runtime = 30
maxTurnaround = Int64(45/5)
top_c = 4

(best_obj, best_sol, iterations) = Heuristic2(data, rt, maxTurnaround, runtime, top_c)

println("Heuristic ran $(iterations) iterations")
println("Best solution:")
println("Objective Value: $(best_obj)")
print_chromosome_table(best_sol)

# start_time = time()
# Best_sols = generate_best_initial_solutions(data, rt, top_k = nr)

# Org_obj = Best_sols[nr].fitness

# evtols_init = Best_sols[nr].evtols

# println("Initial solution:")
# print_chromosome_table(evtols_init)

# # (best_obj, best_sol) = DestructLoop(evtols_init, maxTurnaround, Org_obj, data, rt)
# (best_obj, best_sol) = ConstructLoop(evtols_init, maxTurnaround, Org_obj, data, rt)
# # (best_obj, best_sol) = Swap(evtols_init, maxTurnaround, Org_obj, data, rt)

# println("Best solution:")
# print_chromosome_table(best_sol)

# println("Initial Obj $(Org_obj)")
# println("Best Obj $(best_obj)")