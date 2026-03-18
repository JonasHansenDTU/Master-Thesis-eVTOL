###############################################################################
# UAM / eVTOL optimization model in Julia + JuMP
#
# Reads input from one Excel file with two sheets:
#   1) "Infrastructure"
#   2) "PassengerGroups"

using JuMP
using Gurobi
using XLSX
using DataFrames
using CSV
using MathOptInterface
using Printf
const MOI = MathOptInterface

###############################################################################
# Helpers
###############################################################################

"Normalize column names: lowercase, strip spaces, replace spaces with underscores."
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
    pax   = read_sheet(excel_file, "PassengerGroups")
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

    vb = Dict{Int,Int}()  # base vertiport for each eVTOL
    for r in eachrow(plane)
        n = Int(r[plane_id_col])
        b = Int(r[base_vp_col])
        vb[n] = b
    end

    
    bad_bases = sort([b for b in values(vb) if !(b in V)])
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
    L1                    = 1
    L2                    = bmax*2
    L3                    = ET


    ###########################################################################
    # Node coordinates and parking capacities
    ###########################################################################
    lat = Dict{Int,Float64}()
    lon = Dict{Int,Float64}()
    cap_node = Dict{Int,Int}()

    for r in eachrow(infra)
        j = Int(r[id_col])
        cap_node[j] = Int(r[pads_col])

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
        vb = vb,
        lat = lat, lon = lon,
        dist = dist, fd = fd, fs = fs, c = c, e = e, rt = rt,
        op = op, dp = dp, dt = dt, q = q, so = so, p = p, d = d,
        cap_node = cap_node, cap_flt = cap_flt, cap_u = cap_u,
        bmax = bmax, bmin = bmin, ec = ec, te = te, w = w, ET = ET, L1 = L1, L2 = L2, L3 = L3
    )
end

###############################################################################
# Model builder
###############################################################################

function build_model(excel_file::String; show_progress::Bool = true, display_interval_sec::Int = 5)

    data = load_data(excel_file)

    A = data.A
    V = data.V
    N = data.N
    M = data.M
    M_no0 = data.M_no0
    M_mid = data.M_mid
    M_no_last = data.M_no_last
    T = data.T
    T_no0 = data.T_no0

    vb = data.vb
    fd = data.fd
    fs = data.fs
    c  = data.c
    e  = data.e
    rt = data.rt
    d  = data.d
    op = data.op
    dp = data.dp
    dt = data.dt
    q  = data.q
    so = data.so
    p  = data.p

    cap_node = data.cap_node
    cap_flt  = data.cap_flt
    cap_u    = data.cap_u
    bmax     = data.bmax
    bmin     = data.bmin
    ec       = data.ec
    te       = data.te
    w        = data.w
    L1        = data.L1
    L2        = data.L2
    L3        = data.L3

    ###########################################################################
    # Solver
    ###########################################################################

    model = Model(Gurobi.Optimizer)
    # Show Gurobi MIP progress (incumbent, bound, gap, nodes, time).
    set_optimizer_attribute(model, "OutputFlag", show_progress ? 1 : 0)
    if show_progress
        set_optimizer_attribute(model, "DisplayInterval", max(1, display_interval_sec))
    end

    ###########################################################################
    # Decision variables
    ###########################################################################

    # x[i,j,m,n] = 1 if eVTOL n travels from i to j in operation m
    @variable(model, x[i in V, j in V, m in M, n in N], Bin)

    # s[a,m,n] = 1 if eVTOL n serves passenger group a in operation m
    @variable(model, s[a in A, m in M, n in N], Bin)

    # s[a,n] = 1 if eVTOL n serves passenger group a
    @variable(model, ss[a in A, n in N], Bin)

    # is_p[j,n,t] = 1 if eVTOL n is parking at vertiport/vertistop j at time t
    @variable(model, is_p[j in V, n in N, t in T], Bin)

    # is_o[i,j,m,n,t] = 1 if eVTOL n is traveling from i to j in operation m at time t
    @variable(model, is_o[i in V, j in V, m in M, n in N, t in T], Bin)

    # z[a] = 1 if passenger group a is served by a direct route 1, otherwise 0
    @variable(model, z[a in A], Bin)

    # k[a,i,j,m,n] = 1 if eVTOL n travels from i to j in operation m and serves passenger group a
    @variable(model, k[a in A, i in V, j in V, m in M, n in N], Bin)

    # u[m,n] = battery level of eVTOL n after operation m
    @variable(model, u[m in M, n in N] >= 0)

    # dep[m,n] = departure time of operation m of eVTOL n
    @variable(model, dep[m in M, n in N] >= 0)

    # arr[m,n] = arrival time of operation m of eVTOL n
    @variable(model, arr[m in M, n in N] >= 0)

    ###########################################################################
    # Initialization helpers (operation 0 should not be an actual flown trip)
    ###########################################################################
    for i in V, j in V, n in N
        fix(x[i,j,0,n], 0.0; force=true)
        for a in A
            fix(k[a,i,j,0,n], 0.0; force=true)
        end
    end
    for a in A, n in N
        fix(s[a,0,n], 0.0; force=true)
    end
    
    ###########################################################################
    # Objective (6.1)
    ###########################################################################
    @objective(
        model, Max,
        sum(d[(a,i,j)] * ss[a,n] * (fd[i,j]*(1 - so[a])+ fs[i,j]* so[a]) for a in A, i in V, j in V, n in N) -
        sum(c[(i,j)] * x[i,j,m,n] for i in V, j in V, m in M, n in N) -
        sum(p[a] * (1 - sum(ss[a,n] for n in N)) for a in A)
    )

    ###########################################################################
    # Constraints
    ###########################################################################

    # (6.2) eVTOL leaves base vertiport at time/operation 1
    @constraint(model, [n in N], sum(x[vb[n], j, 1, n] for j in V) == 1)

    # (6.3) eVTOL returns to its base vertiport
    @constraint(model, [n in N],
        sum(x[vb[n], j, m, n] for j in V, m in M_no0) ==
        sum(x[j, vb[n], m, n] for j in V, m in M_no0)
    )

    # (6.4) In each operation m and eVTOL n, at most one trip can be serviced
    @constraint(model, [m in M_no0, n in N],
        sum(x[i,j,m,n] for i in V, j in V) <= 1
    )

    # (6.5) Flow consistency between operation m and m+1
    @constraint(model, [j in V, m in M_mid, n in N],
        sum(x[i,j,m,n] for i in V)  >= sum(x[j,i2,m+1,n] for i2 in V)
    )

    # (6.6) One eVTOL can only contain 2 different passenger groups at the same time
    @constraint(model, [m in M_no0, n in N],
        sum(s[a,m,n] for a in A) <= 2
    )

    # (6.7) Number of flight leg assigned to passenger group a = service indicator + stop indicator 
    @constraint(model, [a in A],
        sum(s[a,m,n] for m in M_no0, n in N) == sum(ss[a,n] for n in N) + z[a]
    )

    # (6.8) Direct connection if z[a] = 0
    @constraint(model, [a in A, m in M, n in N],
        s[a,m,n] <= x[op[a], dp[a], m, n] + z[a] * L1
    )

    # (6.9) Layover path existence upper bound
    @constraint(model, [a in A, m in M, n in N],
        s[a,m,n] <=
        sum(x[op[a], k_node, m, n] for k_node in V) +
        sum(x[k_node, dp[a], m, n] for k_node in V) +
        (1 - z[a]) * L1
    )

    # (6.10) Layover path existence lower bound using m and m+1
    @constraint(model, [a in A, m in M_no_last, n in N],
        s[a,m,n] >=
        sum(x[op[a], k_node, m, n] + x[k_node, dp[a], m+1, n] for k_node in V) - 1 - (1 - z[a]) * L1
    )

    # (6.11) Same as above, using m-1 and m
    @constraint(model, [a in A, m in M_no0, n in N],
        s[a,m,n] >=
        sum(x[op[a], k_node, m-1, n] + x[k_node, dp[a], m, n] for k_node in V) - 1 - (1 - z[a]) * L1
    )

    # (6.12) If passenger group does not allow layovers, then z[a] must be 0
    @constraint(model, [a in A], z[a] <= so[a])

    # (6.13) If passenger group a is not served by eVTOL n, it cannot be served in any operation m
    @constraint(model, [a in A, m in M, n in N],
        s[a,m,n] <= ss[a,n]
    )

    # (6.14) Each passenger group can only be served by one eVTOL
    @constraint(model, [a in A],
        sum(ss[a,n] for n in N) <= 1
    )

    # (6.15) If s[a,m,n] = 1, then at least one arc serving that passenger must exist
    @constraint(model, [a in A, m in M_no0, n in N],
        sum(k[a,i,j,m,n] for i in V, j in V) >= s[a,m,n]
    )

    # (6.16) Direct service linkage
    @constraint(model, [a in A, i in V, j in V, m in M_no0, n in N],
        2 * k[a,i,j,m,n] <= d[(a,i,j)] + x[i,j,m,n] + z[a] * L1
    )

    # (6.17) Layover service linkage
    @constraint(model, [a in A, i in V, k_node in V, j in V, m in M_no_last, n in N],
        k[a,i,k_node,m,n] + k[a,k_node,j,m+1,n] <=
        d[(a,i,j)] + x[i,k_node,m,n] + x[k_node,j,m+1,n] + (1 - z[a]) * L1
    )

    # (6.18) Seat capacity
    @constraint(model, [m in M, n in N],
        sum(s[a,m,n] * q[a] for a in A) <= cap_u
    )

    # (6.19) eVTOL starts with max battery
    @constraint(model, [n in N], u[0,n] == bmax)

    # (6.20) Battery cannot exceed max
    @constraint(model, [m in M, n in N], u[m,n] <= bmax)

    # (6.21) Battery must stay above minimum
    @constraint(model, [m in M, n in N], u[m,n] >= bmin)

    # (6.22) First operation from a vertiport only reflects energy consumption
    @constraint(model, [i in V, j in V, n in N],
        u[1,n] <= u[0,n] - e[(i,j)] * x[i,j,1,n] + (1 - x[i,j,1,n]) * L2
    )

    @constraint(model, [i in V, j in V, n in N],
        u[1,n] >= u[0,n] - e[(i,j)] * x[i,j,1,n] - (1 - x[i,j,1,n]) * L2
    )

    # (6.23) Battery update between operations
    @constraint(model, [i in V, j in V, m in 2:maximum(M), n in N],
        u[m,n] <= u[m-1,n] - e[(i,j)] * x[i,j,m,n] +
                  ec * (arr[m,n] - arr[m-1,n] - rt[(i,j)]) +
                  (1 - x[i,j,m,n]) * L2
    )

    @constraint(model, [i in V, j in V, m in 2:maximum(M), n in N],
        u[m,n] >= u[m-1,n] - e[(i,j)] * x[i,j,m,n] +
                  ec * (arr[m,n] - arr[m-1,n] - rt[(i,j)]) -
                  (1 - x[i,j,m,n]) * L2
    )

    # (6.24) Operation 0 starts at time 0
    @constraint(model, [n in N], arr[0,n] == 0)

    # (6.25) Arrival time lower bound
    @constraint(model, [m in M_no0, n in N],
    arr[m,n] >= arr[m-1,n] + sum((te + rt[(i,j)]) * x[i,j,m,n] for i in V, j in V)
    )

    # (6.26) Departure time = arrival time - travel time
    @constraint(model, [m in M, n in N],
    dep[m,n] == arr[m,n] - sum(rt[(i,j)] * x[i,j,m,n] for i in V, j in V)
    )

    # (6.27) Minimum layover time
    @constraint(model, [a in A, n in N, m in M_no_last],
        dep[m+1,n] <= arr[m,n] + te + (2 - s[a,m,n] - s[a,m+1,n]) * L3
    )

    # (6.28) eVTOL must arrive before passenger group arrives for boarding
    @constraint(model, [a in A, m in M_no0, n in N],
        arr[m-1,n] <= sum(d[(a,i,j)] * dt[a] for i in V, j in V) +
                      (1 - (s[a,m,n] - s[a,m-1,n])) * L3
    )

    # (6.29) Earliest arrival time at destination
    @constraint(model, [a in A, i in V, j in V, m in M_no0, n in N],
        d[(a,i,j)] * dt[a] - (1 - (s[a,m,n] - s[a,m-1,n])) * L3 <=
        arr[m,n] - sum(rt[(i,k_node)] * x[i,k_node,m,n] for k_node in V)
    )

    # (6.30) Maximum waiting time
    @constraint(model, [a in A, m in M_no0, n in N],
        arr[m,n]
        - sum(rt[(i,j)] * x[i,j,m,n] for i in V, j in V)
        - sum(d[(a,i,j)] * dt[a] * s[a,m,n] for i in V, j in V)
        - (1 - (s[a,m,n] - s[a,m-1,n])) * L3 <= w
    )

    # (6.31) eVTOL is either parked or flying at each time t
    @constraint(model, [n in N, t in T],
        sum(is_p[j,n,t] for j in V) +
        sum(is_o[i,j,m,n,t] for i in V, j in V, m in M) == 1
    )

    # (6.32) Travel time occupancy relation
    @constraint(model, [i in V, j in V, m in M, n in N],
        rt[(i,j)] * x[i,j,m,n] == sum(is_o[i,j,m,n,t] for t in T)
    )

    # (6.33) Departure time bound from occupancy
    @constraint(model, [i in V, j in V, m in M_no0, n in N, t in T],
        dep[m,n] <= t + L3 * (1 - is_o[i,j,m,n,t]) - 1
    )

    # (6.34) Arrival time bound from occupancy
    @constraint(model, [i in V, j in V, m in M_no0, n in N, t in T],
        arr[m,n] >= t - L3 * (1 - is_o[i,j,m,n,t])
    )

    # (6.35) Parking state propagation
    @constraint(model, [j in V, n in N, t in T_no0],
        is_p[j,n,t] <= is_p[j,n,t-1] +
                       sum(is_o[i,j,m,n,t-1] for i in V, m in M)
    )

    # (6.36) Parking capacity at vertiports
    @constraint(model, [j in V, t in T],
        sum(is_p[j,n,t] for n in N) <= cap_node[j]
    )

    # (6.37) Air corridor capacity
    @constraint(model, [i in V, j in V, t in T],
        sum(is_o[i,j,m,n,t] for m in M, n in N) <= cap_flt
    )

    # (6.38) Initial parking at base vertiport
    @constraint(model, [n in N],
        is_p[vb[n], n, 0] == 1
    )

    # (6.39) 
    @constraint(model, [i in V, m in M_no0, n in N], x[i,i,m,n] <= 0)

    return model, data
end

###############################################################################
# Solve + simple reporting
###############################################################################

function solve_instance(excel_file::String; show_progress::Bool = true, display_interval_sec::Int = 5)
    timings = Dict{String,Float64}()

    t_build = @elapsed model, data = build_model(
        excel_file;
        show_progress = show_progress,
        display_interval_sec = display_interval_sec,
    )
    timings["Build model (incl. data load)"] = t_build

    t_opt = @elapsed optimize!(model)
    timings["Optimization"] = t_opt

    term = termination_status(model)
    primal = primal_status(model)

    println("Termination status: ", term)
    println("Primal status:      ", primal)

    if term == MOI.OPTIMAL || term == MOI.TIME_LIMIT || term == MOI.FEASIBLE_POINT
        println("Objective value:    ", objective_value(model))
    end

    return model, data, timings
end

function print_timing_summary(timings::Dict{String,Float64})
    println("\n")
    println("=" ^ 80)
    println("RUNTIME SUMMARY")
    println("=" ^ 80)
    println(lpad("Program Part", 38) * " | " * lpad("Seconds", 10))
    println("-" ^ 53)

    measured_parts = [
        "Build model (incl. data load)",
        "Optimization",
        "Snapshot export",
        "Pretty printing"
    ]

    subtotal = 0.0
    for part in measured_parts
        if haskey(timings, part)
            v = timings[part]
            subtotal += v
            println(lpad(part, 38) * " | " * lpad(@sprintf("%.3f", v), 10))
        end
    end

    if haskey(timings, "Total script")
        total = timings["Total script"]
        overhead = total - subtotal
        println("-" ^ 53)
        println(lpad("Measured subtotal", 38) * " | " * lpad(@sprintf("%.3f", subtotal), 10))
        println(lpad("Unaccounted overhead", 38) * " | " * lpad(@sprintf("%.3f", overhead), 10))
        println(lpad("Total script", 38) * " | " * lpad(@sprintf("%.3f", total), 10))
    end
end

function export_solution_snapshots(model::Model, data; out_csv::String = joinpath(@__DIR__, "solution_snapshots.csv"))
    if !has_values(model)
        println("No primal solution available; skipping snapshot export.")
        return nothing
    end

    V = data.V
    A = data.A
    N = data.N
    M = data.M
    T = data.T
    q = data.q
    lat = data.lat
    lon = data.lon

    # Precompute passenger groups served by each eVTOL (global assignment).
    served_groups_by_evtol = Dict{Int,Vector{Int}}()
    for n in N
        served_groups_by_evtol[n] = [a for a in A if value(model[:ss][a, n]) > 0.5]
    end

    rows = NamedTuple[]

    for n in N, t in T
        # Best parked state at (n, t)
        best_p = 0.0
        park_node = missing
        for j in V
            pv = value(model[:is_p][j, n, t])
            if pv > best_p
                best_p = pv
                park_node = j
            end
        end

        # Best flying state at (n, t)
        best_o = 0.0
        fly_i = missing
        fly_j = missing
        fly_m = missing
        for i in V, j in V, m in M
            ov = value(model[:is_o][i, j, m, n, t])
            if ov > best_o
                best_o = ov
                fly_i = i
                fly_j = j
                fly_m = m
            end
        end

        state = best_o > 0.5 ? "flying" : (best_p > 0.5 ? "parked" : "inactive")

        node_from = state == "flying" ? fly_i : (state == "parked" ? park_node : missing)
        node_to = state == "flying" ? fly_j : (state == "parked" ? park_node : missing)
        op = state == "flying" ? fly_m : missing

        x_from = (node_from === missing) ? missing : lon[node_from]
        y_from = (node_from === missing) ? missing : lat[node_from]
        x_to = (node_to === missing) ? missing : lon[node_to]
        y_to = (node_to === missing) ? missing : lat[node_to]

        x = (x_from === missing || x_to === missing) ? missing : (x_from + x_to) / 2
        y = (y_from === missing || y_to === missing) ? missing : (y_from + y_to) / 2

        # Battery estimate at time t: use battery level after latest arrived operation.
        latest_completed_m = minimum(M)
        for m in sort(M)
            if value(model[:arr][m, n]) <= t + 1e-6
                latest_completed_m = m
            end
        end
        battery_level = value(model[:u][latest_completed_m, n])

        # Passenger load at time t: if flying in op m, use groups served in that op.
        onboard_groups = Int[]
        if state == "flying" && fly_m !== missing
            onboard_groups = [a for a in A if value(model[:s][a, fly_m, n]) > 0.5]
        end
        onboard_passenger_count = sum((q[a] for a in onboard_groups); init=0)
        onboard_group_sizes = isempty(onboard_groups) ? "" : join(["$(a):$(q[a])" for a in onboard_groups], ";")

        served_groups_evtol = served_groups_by_evtol[n]

        push!(rows, (
            time = t,
            evtol_id = n,
            state = state,
            node_from = node_from,
            node_to = node_to,
            op = op,
            is_p = best_p,
            is_o = best_o,
            battery_level = battery_level,
            battery_after_op = latest_completed_m,
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

    snapshots = DataFrame(rows)
    
    # Append infrastructure nodes as vertiport markers
    for j in V
        push!(snapshots, (
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
    
    CSV.write(out_csv, snapshots)
    println("Snapshot export written: ", out_csv, " (rows=", nrow(snapshots), ")")
    return snapshots
end

###############################################################################
# Usage
###############################################################################

excel_file = joinpath(@__DIR__, "inputData.xlsx")
println("Using Excel file: ", excel_file)
total_start = time()
model, data, timings = solve_instance(excel_file)

t_export = @elapsed export_solution_snapshots(model, data)
timings["Snapshot export"] = t_export

###############################################################################
# Pretty result printing - ALL VARIABLES
###############################################################################

function print_results_pretty(model::Model, data)
    divider(w=80) = println("=" ^ w)
    section(title) = begin
        println("\n")
        divider()
        println(title)
        divider()
    end

    A = data.A
    V = data.V
    N = data.N
    M = data.M
    M_no0 = data.M_no0
    M_mid = data.M_mid
    T = data.T
    T_no0 = data.T_no0
    vb = data.vb
    cap_node = data.cap_node
    cap_flt = data.cap_flt
    cap_u = data.cap_u
    bmax = data.bmax
    bmin = data.bmin
    ec = data.ec
    te = data.te
    w = data.w
    ET = data.ET
    L1 = data.L1
    L2 = data.L2
    L3 = data.L3
    fd = data.fd
    fs = data.fs
    c = data.c
    e = data.e
    rt = data.rt
    rt_int = hasproperty(data, :rt_int) ? data.rt_int : Dict((i, j) => Int(ceil(rt[(i, j)])) for i in V, j in V)
    op = data.op
    dp = data.dp
    dt = data.dt
    q = data.q
    so = data.so
    p = data.p
    d = data.d
    dist = data.dist

    section("SETS")
    println("Vertices (V): ", V)
    println("Passenger groups (A): ", A)
    println("eVTOLs (N): ", N)
    println("Operations (M): ", M)
    println("Operations (M_no0): ", M_no0)
    println("Operations (M_mid): ", M_mid)
    println("Time periods (T): 0 to ", maximum(T))
    println("Time periods (T_no0): 1 to ", maximum(T_no0))
    println("Base vertiports (vb): ", vb)

    section("SYSTEM PARAMETERS")
    println(lpad("Battery capacity (bmax)", 35) * ": " * @sprintf("%.2f", bmax) * " %")
    println(lpad("Battery minimum (bmin)", 35) * ": " * @sprintf("%.2f", bmin) * " %")
    println(lpad("Charging rate (ec)", 35) * ": " * @sprintf("%.2f", ec) * " %/unit")
    println(lpad("Turnaround time (te)", 35) * ": " * @sprintf("%.2f", te) * " min")
    println(lpad("Maximum waiting time (w)", 35) * ": " * @sprintf("%.2f", w) * " min")
    println(lpad("End time (ET)", 35) * ": " * @sprintf("%.2f", ET))
    println(lpad("eVTOL seat capacity (cap_u)", 35) * ": " * string(Int(cap_u)))
    println(lpad("Air corridor capacity (cap_flt)", 35) * ": " * string(Int(cap_flt)))

    section("INFRASTRUCTURE NODES & PARKING CAPACITY")
    println(lpad("Node", 6) * " | " * lpad("Type", 11) * " | " * lpad("Capacity", 10))
    println("-" ^ 32)

    section("DISTANCES (dist[i,j])")
    println(lpad("From", 6) * " | " * lpad("To", 6) * " | " * lpad("Distance (km)", 15))
    println("-" ^ 32)
    for i in V, j in V
        if i <= j
            println(lpad(string(i), 6) * " | " * lpad(string(j), 6) * " | " * lpad(@sprintf("%.2f", dist[(i,j)]), 15))
        end
    end

    section("ARC PARAMETERS: fd, fs, c, e, rt, rt_int")
    println(lpad("i->j", 5) * " | " * lpad("Fare Direct", 14) * " | " * lpad("Fare Stop", 12) *
            " | " * lpad("Cost", 7) * " | " * lpad("Battery", 9) * " | " * lpad("Time", 7) * " | " * lpad("Time_int", 9))
    println("-" ^ 85)
    for i in V, j in V
        if i != j
            println(lpad("$i->$j", 5) * " | " * lpad(@sprintf("%.2f", fd[(i,j)]), 14) * " | " *
                    lpad(@sprintf("%.2f", fs[(i,j)]), 12) * " | " * lpad(@sprintf("%.2f", c[(i,j)]), 7) *
                    " | " * lpad(@sprintf("%.2f", e[(i,j)]), 9) * " | " * lpad(@sprintf("%.2f", rt[(i,j)]), 7) *
                    " | " * lpad(string(rt_int[(i,j)]), 9))
        end
    end

    section("PASSENGER GROUPS: op, dp, dt, q, so, p")
    println(lpad("G", 4) * " | " * lpad("Origin", 8) * " | " * lpad("Dest", 6) * " | " *
            lpad("Time", 8) * " | " * lpad("Pax", 5) * " | " * lpad("Stopover", 10) * " | " * lpad("Penalty", 9))
    println("-" ^ 68)
    for a in A
        println(lpad(string(a), 4) * " | " * lpad(string(op[a]), 8) * " | " * lpad(string(dp[a]), 6) * " | " *
                lpad(@sprintf("%.0f", dt[a]), 8) * " | " * lpad(string(q[a]), 5) * " | " *
                lpad(so[a] > 0.5 ? "Yes" : "No", 10) * " | " * lpad(@sprintf("%.2f", p[a]), 9))
    end

    section("DESTINATION INDICATOR: d[a,i,j]")
    println(lpad("Group", 6) * " | " * lpad("Origin", 8) * " | " * lpad("Dest", 6) * " | " * lpad("d[a,i,j]", 10))
    println("-" ^ 35)
    for a in A, i in V, j in V
        if d[(a,i,j)] > 0
            println(lpad(string(a), 6) * " | " * lpad(string(i), 8) * " | " * lpad(string(j), 6) * " | " * lpad(string(d[(a,i,j)]), 10))
        end
    end

    if !has_values(model)
        println("\n")
        divider()
        println("No solution available - skipping decision variables")
        divider()
        return
    end

    section("DECISION VAR: x[i,j,m,n] (Flight Arcs)")
    count = 0
    for n in N, m in M, i in V, j in V
        if value(model[:x][i,j,m,n]) > 0.5
            if count == 0
                println(lpad("eVTOL", 8) * " | " * lpad("Op", 5) * " | " * lpad("From", 6) * " | " * lpad("To", 6) * " | " * lpad("Value", 8))
                println("-" ^ 42)
            end
            println(lpad(string(n), 8) * " | " * lpad(string(m), 5) * " | " * lpad(string(i), 6) * " | " * lpad(string(j), 6) * " | " * lpad(@sprintf("%.2f", value(model[:x][i,j,m,n])), 8))
            count += 1
        end
    end
    if count == 0
        println("(all zeros)")
    end

    section("DECISION VAR: s[a,m,n] (Service per Operation)")
    println(lpad("Group", 6) * " | " * lpad("Op", 5) * " | " * lpad("eVTOL", 8) * " | " * lpad("Value", 8))
    println("-" ^ 35)
    for a in A, m in M_no0, n in N
        sval = value(model[:s][a,m,n])
        if sval > 0.01
            println(lpad(string(a), 6) * " | " * lpad(string(m), 5) * " | " * lpad(string(n), 8) * " | " * lpad(@sprintf("%.2f", sval), 8))
        end
    end

    section("DECISION VAR: ss[a,n] (Service by eVTOL)")
    println(lpad("Group", 6) * " | " * lpad("eVTOL", 8) * " | " * lpad("Value", 8))
    println("-" ^ 28)
    for a in A, n in N
        sval = value(model[:ss][a,n])
        if sval > 0.01
            println(lpad(string(a), 6) * " | " * lpad(string(n), 8) * " | " * lpad(@sprintf("%.2f", sval), 8))
        end
    end

    section("DECISION VAR: z[a] (Direct=0, Layover=1)")
    println(lpad("Group", 6) * " | " * lpad("Value", 8) * " | " * lpad("Type", 12))
    println("-" ^ 32)
    for a in A
        zval = value(model[:z][a])
        ztype = zval > 0.5 ? "Layover" : "Direct"
        println(lpad(string(a), 6) * " | " * lpad(@sprintf("%.2f", zval), 8) * " | " * lpad(ztype, 12))
    end

    section("DECISION VAR: k[a,i,j,m,n] (Service Arc)")
    count = 0
    for a in A, i in V, j in V, m in M_no0, n in N
        kval = value(model[:k][a,i,j,m,n])
        if kval > 0.5
            if count == 0
                println(lpad("Group", 6) * " | " * lpad("i->j", 5) * " | " * lpad("Op", 5) * " | " * lpad("eVTOL", 8) * " | " * lpad("Value", 8))
                println("-" ^ 42)
            end
            println(lpad(string(a), 6) * " | " * lpad("$i->$j", 5) * " | " * lpad(string(m), 5) * " | " * lpad(string(n), 8) * " | " * lpad(@sprintf("%.2f", kval), 8))
            count += 1
        end
    end
    if count == 0
        println("(all zeros)")
    end

    section("DECISION VAR: u[m,n] (Battery Level %)")
    println(lpad("eVTOL", 8) * " | " * lpad("Op", 5) * " | " * lpad("Battery %", 12))
    println("-" ^ 32)
    for n in N, m in M
        u_val = value(model[:u][m,n])
        println(lpad(string(n), 8) * " | " * lpad(string(m), 5) * " | " * lpad(@sprintf("%.2f", u_val), 12))
    end

    section("DECISION VAR: dep[m,n] and arr[m,n] (Times)")
    println(lpad("eVTOL", 8) * " | " * lpad("Op", 5) * " | " * lpad("Dep (min)", 12) * " | " * lpad("Arr (min)", 12))
    println("-" ^ 48)
    for n in N, m in M
        dep_val = value(model[:dep][m,n])
        arr_val = value(model[:arr][m,n])
        println(lpad(string(n), 8) * " | " * lpad(string(m), 5) * " | " * lpad(@sprintf("%.2f", dep_val), 12) * " | " * lpad(@sprintf("%.2f", arr_val), 12))
    end

    section("DECISION VAR: is_p[j,n,t] and is_o[i,j,m,n,t] (Occupancy States)")
    println(lpad("eVTOL", 8) * " | " * lpad("Time", 6) * " | " * lpad("Park Node", 10) * " | " * lpad("Travel Arc", 10) * " | " * lpad("Op", 5) * " | " * lpad("is_p", 8) * " | " * lpad("is_o", 8))
    println("-" ^ 73)

    shown = 0
    p_count = 0
    o_count = 0

    for n in N, t in T
        p_node = "-"
        check = 0.5
        p_val = 0.0
        for j in V
            pv = value(model[:is_p][j,n,t])
            if pv > check
                p_val = pv
                p_node = string(j)
            end
        end

        o_arc = "-"
        o_op = "-"
        o_val = 0.0
        for i in V, j in V, m in M
            ov = value(model[:is_o][i,j,m,n,t])
            if ov > check
                o_val = ov
                o_arc = "$i->$j"
                o_op = string(m)
            end
        end

        if p_val > 0.01 || o_val > 0.01
            println(lpad(string(n), 8) * " | " * lpad(string(t), 6) * " | " * lpad(p_node, 10) * " | " * lpad(o_arc, 10) * " | " * lpad(o_op, 5) * " | " * lpad(@sprintf("%.2f", p_val), 8) * " | " * lpad(@sprintf("%.2f", o_val), 8))
            shown += 1
            if p_val > 0.5
                p_count += 1
            end
            if o_val > 0.5
                o_count += 1
            end
        end
    end

    if shown == 0
        println("(all or mostly zeros)")
    end
    println("Summary: rows shown = $(shown), active is_p = $(p_count), active is_o = $(o_count)")

    divider()
    println()
end

t_pretty = @elapsed Base.invokelatest(print_results_pretty, model, data)
timings["Pretty printing"] = t_pretty
timings["Total script"] = time() - total_start

print_timing_summary(timings)