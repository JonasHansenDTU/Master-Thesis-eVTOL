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
using MathOptInterface
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
Build all sets and parameters from Excel + hard-coded system parameters.
Expected sheets:
  - Infrastructure
  - PassengerGroups
"""
function load_data(excel_file::String)

    ###########################################################################
    # Read sheets
    ###########################################################################
    infra = read_sheet(excel_file, "Infrastructure")
    pax   = read_sheet(excel_file, "PassengerGroups")

    ###########################################################################
    # Infrastructure columns
    ###########################################################################
    type_col   = find_col(infra, [:type])
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

    VP = sort(Int.(infra[lowercase.(String.(infra[!, type_col])) .== "vertiport", id_col]))
    VS = sort(Int.(infra[lowercase.(String.(infra[!, type_col])) .== "vertistop", id_col]))

    N = 1:4
    vb = Dict(1 => 1, 2 => 1, 3 => 2, 4 => 2)

    M = 0:5
    M_no0 = 1:maximum(M)
    M_mid = 1:(maximum(M)-1)

    T = 0:120
    T_no0 = 1:maximum(T)

    # Passenger groups
    A = sort(Int.(pax[!, group_col]))

    ###########################################################################
    # System parameters from Table 3
    ###########################################################################
    fare_direct_per_km    = 3.0      # fd_{i,j} = $3/km
    fare_stopover_factor  = 0.75     # fs_{i,j} = 75% of fd_{i,j}
    p_penalty             = 10.0     # p_a = $10
    operating_cost_per_km = 1.0      # c_{i,j} = $1/km
    battery_per_km        = 1.0      # e_{i,j} = 1%/km
    time_per_km           = 0.6      # rt_{i,j} = 0.6 min/km
    cap_flt               = 1        # capacity of air corridor
    cap_u                 = 4        # seat capacity of eVTOL
    bmax                  = 100.0    # 100%
    bmin                  = 10.0     # 10%
    ec                    = 5.0      # 5% battery charged per unit time
    te                    = 10.0     # minimum turnaround time at vertiport
    w                     = 10.0     # maximum waiting time
    ET                    = 500      # end time
    L                     = 10_000.0 # default Big-M

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
    rt   = Dict{Tuple{Int,Int},Float64}()
    rt_int = Dict{Tuple{Int,Int},Int}()

    for i in V, j in V
        dij = haversine_km(lat[i], lon[i], lat[j], lon[j])
        dist[(i,j)] = dij
        fd[(i,j)] = fare_direct_per_km * dij
        fs[(i,j)] = fare_stopover_factor * fd[(i,j)]
        c[(i,j)]  = operating_cost_per_km * dij
        e[(i,j)]  = battery_per_km * dij
        rt[(i,j)] = time_per_km * dij
        rt_int[(i,j)] = Int(ceil(rt[(i,j)]))   # discretized travel time for time-indexed constraints
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
        V = V, VP = VP, VS = VS, A = A, N = collect(N),
        M = collect(M), M_no0 = collect(M_no0), M_mid = collect(M_mid),
        T = collect(T), T_no0 = collect(T_no0),
        vb = vb,
        dist = dist, fd = fd, fs = fs, c = c, e = e, rt = rt, rt_int = rt_int,
        op = op, dp = dp, dt = dt, q = q, so = so, p = p, d = d,
        cap_node = cap_node, cap_flt = cap_flt, cap_u = cap_u,
        bmax = bmax, bmin = bmin, ec = ec, te = te, w = w, ET = ET, L = L
    )
end

###############################################################################
# Model builder
###############################################################################

function build_model(excel_file::String)

    data = load_data(excel_file)

    A = data.A
    V = data.V
    VP = data.VP
    VS = data.VS
    N = data.N
    M = data.M
    M_no0 = data.M_no0
    M_mid = data.M_mid
    T = data.T
    T_no0 = data.T_no0

    vb = data.vb
    fd = data.fd
    fs = data.fs
    c  = data.c
    e  = data.e
    rt = data.rt
    rt_int = data.rt_int
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
    L        = data.L

    ###########################################################################
    # Solver
    ###########################################################################

    model = Model(Gurobi.Optimizer)
    set_silent(model)

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
        sum(d[(a,i,j)] * ss[a,n] * (fd[i,j]*(1 - so[a]) + fs[i,j]*so[a]) for a in A, i in V, j in V, n in N) -
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
        sum(x[i,j,m,n] for i in V) - sum(x[j,i2,m+1,n] for i2 in V) >= 0
    )

    # (6.6) No self-travel i -> i
    @constraint(model, [i in V, m in M_no0, n in N], x[i,i,m,n] <= 0)

    # (6.7) Once an eVTOL stops flying, it cannot restart later
    @constraint(model, [m in M_mid, n in N],
        sum(x[i,j,m,n] for i in V, j in V) >=
        sum(x[i,j,m+1,n] for i in V, j in V)
    )

    # (6.8) One eVTOL can only contain 2 different passenger groups at the same time
    @constraint(model, [m in M_no0, n in N],
        sum(s[a,m,n] for a in A) <= 2
    )

    # (6.9) One passenger group can only be part of 2 trips per eVTOL
    @constraint(model, [a in A, n in N],
        sum(s[a,m,n] for m in M_no0) <= 2
    )

    # (6.10)
    @constraint(model, [a in A],
        sum(s[a,m,n] for m in M_no0, n in N) <= 2*z[a] + 1
    )

    # (6.11)
    @constraint(model, [a in A],
        sum(s[a,m,n] for m in M_no0, n in N) >= 2*z[a]
    )

    # (6.12) Direct connection if z[a] = 0
    @constraint(model, [a in A, m in M, n in N],
        s[a,m,n] <= x[op[a], dp[a], m, n] + z[a] * L
    )

    # (6.13) Layover path existence upper bound
    @constraint(model, [a in A, m in M, n in N],
        s[a,m,n] <=
        sum(x[op[a], k_node, m, n] for k_node in V) +
        sum(x[k_node, dp[a], m, n] for k_node in V) +
        (1 - z[a]) * L
    )

    # (6.14) Layover path existence lower bound using m and m+1
    @constraint(model, [a in A, m in M_mid, n in N],
        s[a,m,n] >=
        sum(x[op[a], k_node, m, n] + x[k_node, dp[a], m+1, n] for k_node in V) - 1 - (1 - z[a]) * L
    )

    # (6.15) Same as above, using m-1 and m
    @constraint(model, [a in A, m in M_no0, n in N],
        s[a,m,n] >=
        sum(x[op[a], k_node, m-1, n] + x[k_node, dp[a], m, n] for k_node in V) - 1 - (1 - z[a]) * L
    )

    # (6.16) If passenger group does not allow layovers, then z[a] must be 0
    @constraint(model, [a in A], z[a] <= so[a])

    # (6.17) If passenger group a is not served by eVTOL n, it cannot be served in any operation m
    @constraint(model, [a in A, m in M, n in N],
        s[a,m,n] <= ss[a,n]
    )

    # (6.18) If passenger group a is served by eVTOL n, then it must be served in some operation
    @constraint(model, [a in A, n in N],
        sum(s[a,m,n] for m in M) >= ss[a,n]
    )

    # (6.19) Each passenger group can only be served by one eVTOL
    @constraint(model, [a in A],
        sum(ss[a,n] for n in N) <= 1
    )

    # (6.20) If s[a,m,n] = 1, then at least one arc serving that passenger must exist
    @constraint(model, [a in A, m in M_no0, n in N],
        sum(k[a,i,j,m,n] for i in V, j in V) >= s[a,m,n]
    )

    # (6.21) Direct service linkage
    @constraint(model, [a in A, i in V, j in V, m in M_no0, n in N],
        2 * k[a,i,j,m,n] <= d[(a,i,j)] + x[i,j,m,n] + z[a] * L
    )

    # (6.22) Layover service linkage
    @constraint(model, [a in A, i in V, k_node in V, j in V, m in M_mid, n in N],
        k[a,i,k_node,m,n] + k[a,k_node,j,m+1,n] <=
        d[(a,i,j)] + x[i,k_node,m,n] + x[k_node,j,m+1,n] + (1 - z[a]) * L
    )

    # (6.23) Seat capacity
    @constraint(model, [m in M, n in N],
        sum(s[a,m,n] * q[a] for a in A) <= cap_u
    )

    # (6.24) eVTOL starts with max battery
    @constraint(model, [n in N], u[0,n] == bmax)

    # (6.25) Battery cannot exceed max
    @constraint(model, [m in M, n in N], u[m,n] <= bmax)

    # (6.26) Battery must stay above minimum
    @constraint(model, [m in M, n in N], u[m,n] >= bmin)

    # (6.27) Battery update when departing from a vertistop
    @constraint(model, [i in VS, j in V, m in M_no0, n in N],
        u[m,n] <= u[m-1,n] - e[(i,j)] * x[i,j,m,n] + (1 - x[i,j,m,n]) * L
    )

    # (6.28) First operation from a vertiport only reflects energy consumption
    @constraint(model, [i in VP, j in V, n in N],
        u[1,n] <= u[0,n] - e[(i,j)] * x[i,j,1,n] + (1 - x[i,j,1,n]) * L
    )

    # (6.29) Battery update with charging opportunity between operations
    @constraint(model, [i in VP, j in V, m in 2:maximum(M), n in N],
        u[m,n] <= u[m-1,n] - e[(i,j)] * x[i,j,m,n] +
                  ec * (arr[m,n] - arr[m-1,n] - rt[(i,j)]) +
                  (1 - x[i,j,m,n]) * L
    )

    # (6.30) Operation 0 starts at time 0
    @constraint(model, [n in N], arr[0,n] == 0)

    # (6.31) Arrival time lower bound
    @constraint(model, [i in V, j in V, m in M_no0, n in N],
        arr[m,n] >= arr[m-1,n] + (te + rt[(i,j)]) * x[i,j,m,n]
    )

    # (6.32) Departure time = arrival time - travel time
    @constraint(model, [m in M, n in N],
        dep[m,n] == arr[m,n] - sum(rt[(i,j)] * x[i,j,m,n] for i in V, j in V)
    )

    # (6.33) Minimum layover time
    @constraint(model, [a in A, n in N, m in M_mid],
        dep[m+1,n] <= arr[m,n] + te + (2 - s[a,m,n] - s[a,m+1,n]) * L
    )

    # (6.34) eVTOL must arrive before passenger group arrives for boarding
    @constraint(model, [a in A, m in M_no0, n in N],
        arr[m-1,n] <= sum(d[(a,i,j)] * dt[a] for i in V, j in V) +
                      (1 - (s[a,m,n] - s[a,m-1,n])) * L
    )

    # (6.35) Earliest arrival time at destination
    @constraint(model, [a in A, i in V, j in V, m in M_no0, n in N],
        d[(a,i,j)] * dt[a] - (1 - (s[a,m,n] - s[a,m-1,n])) * L <=
        arr[m,n] - sum(rt[(i,k_node)] * x[i,k_node,m,n] for k_node in V)
    )

    # (6.36) Maximum waiting time
    @constraint(model, [a in A, m in M_no0, n in N],
        arr[m,n]
        - sum(rt[(i,j)] * x[i,j,m,n] for i in V, j in V)
        - sum(d[(a,i,j)] * dt[a] * s[a,m,n] for i in V, j in V)
        - (1 - (s[a,m,n] - s[a,m-1,n])) * L <= w
    )

    # (6.37) eVTOL is either parked or flying at each time t
    @constraint(model, [n in N, t in T],
        sum(is_p[j,n,t] for j in V) +
        sum(is_o[i,j,m,n,t] for i in V, j in V, m in M) == 1
    )

    # (6.38) Travel time occupancy relation
    @constraint(model, [i in V, j in V, m in M, n in N],
        rt_int[(i,j)] * x[i,j,m,n] == sum(is_o[i,j,m,n,t] for t in T)
    )

    # (6.39) Departure time bound from occupancy
    @constraint(model, [i in V, j in V, m in M_no0, n in N, t in T],
        dep[m,n] <= t + L * (1 - is_o[i,j,m,n,t]) - 1
    )

    # (6.40) Arrival time bound from occupancy
    @constraint(model, [i in V, j in V, m in M_no0, n in N, t in T],
        arr[m,n] >= t - L * (1 - is_o[i,j,m,n,t])
    )

    # (6.41) Parking state propagation
    @constraint(model, [j in V, n in N, t in T_no0],
        is_p[j,n,t] <= is_p[j,n,t-1] +
                       sum(is_o[i,j,m,n,t-1] for i in V, m in M)
    )

    # (6.42) Parking capacity at vertiports
    @constraint(model, [j in VP, t in T],
        sum(is_p[j,n,t] for n in N) <= cap_node[j]
    )

    # (6.43) Parking capacity at vertistops
    @constraint(model, [j in VS, t in T],
        sum(is_p[j,n,t] for n in N) <= cap_node[j]
    )

    # (6.44) Air corridor capacity
    @constraint(model, [i in V, j in V, t in T],
        sum(is_o[i,j,m,n,t] for m in M, n in N) <= cap_flt
    )

    # (6.45) Initial parking at base vertiport
    @constraint(model, [n in N],
        is_p[vb[n], n, 0] == 1
    )

    return model, data
end

###############################################################################
# Solve + simple reporting
###############################################################################

function solve_instance(excel_file::String)
    model, data = build_model(excel_file)
    optimize!(model)

    term = termination_status(model)
    primal = primal_status(model)

    println("Termination status: ", term)
    println("Primal status:      ", primal)

    if term == MOI.OPTIMAL || term == MOI.TIME_LIMIT || term == MOI.FEASIBLE_POINT
        println("Objective value:    ", objective_value(model))
    end

    return model, data
end

###############################################################################
# Usage
###############################################################################

excel_file = joinpath(@__DIR__, "inputData.xlsx")
println("Using Excel file: ", excel_file)
model, data = solve_instance(excel_file)

###############################################################################
# Example result extraction
###############################################################################

A = data.A
V = data.V
N = data.N
M = data.M
vb = data.vb

println("\nSelected flight arcs x[i,j,m,n] = 1:")
for n in N, m in M, i in V, j in V
    if value(model[:x][i,j,m,n]) > 0.5
        println("  eVTOL $n, op $m: $i -> $j")
    end
end

println("\nPassenger assignment ss[a,n] = 1:")
for a in A, n in N
    if value(model[:ss][a,n]) > 0.5
        println("  passenger group $a served by eVTOL $n")
    end
end

println("\nDirect/layover indicator z[a]:")
for a in A
    println("  z[$a] = ", value(model[:z][a]))
end

println("\nBattery levels u[m,n]:")
for n in N, m in M
    println("  u[$m,$n] = ", value(model[:u][m,n]))
end

println("\nArrival/departure times:")
for n in N, m in M
    println("  eVTOL $n, op $m: dep = ", value(model[:dep][m,n]),
            ", arr = ", value(model[:arr][m,n]))
end