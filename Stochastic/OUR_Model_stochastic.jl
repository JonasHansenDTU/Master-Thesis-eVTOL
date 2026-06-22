###############################################################################
# OUR_Model_stochastic.jl
#
# Two-stage stochastic extension of the deterministic eVTOL routing model.
#
# FIRST-STAGE (scenario-independent, decided before uncertainty is revealed):
#   y[n]      fleet activation
#   alpha[a]  passenger group acceptance (eVTOL-independent)
#
# SECOND-STAGE (scenario-dependent, decided after scenario is realized):
#   z[sc,a]             stopover indicator (per scenario)
#   x[sc,i,j,m,n]       routing
#   s[sc,a,m,n]          service per operation
#   k[sc,a,i,j,m,n]      service arc
#   u[sc,m,n]            battery level
#   dep[sc,m,n]          departure time
#   arr[sc,m,n]          arrival time
#   charge[sc,m,n]       battery charged during turnaround
#   over_bmid[sc,m,n]    battery above bmid (penalty term)
#   is_p[sc,j,n,t]       parking occupancy
#   is_o[sc,i,j,m,n,t]  flying occupancy
#
# Objective: maximize expected profit across all scenarios.
###############################################################################

using JuMP
using Gurobi
using XLSX
using DataFrames
using CSV
using MathOptInterface
using Printf
using HTTP
using JSON3
using Random
const MOI = MathOptInterface

# Load scenario generation helpers
include("scenario_generation.jl")

###############################################################################
# Helpers  (identical to deterministic model)
###############################################################################

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

function read_sheet(path::String, sheet_name::String)
    df = DataFrame(XLSX.readtable(path, sheet_name, infer_eltypes=true))
    rename!(df, Dict(n => normalize_name(n) for n in names(df)))
    return df
end

function read_sheet_any(path::String, candidates::Vector{String})
    for sheet_name in candidates
        try
            return read_sheet(path, sheet_name)
        catch end
    end
    error("Could not find any of these sheets: $(candidates)")
end

function find_col(df::DataFrame, candidates::Vector{Symbol})
    name_map = Dict(Symbol(String(n)) => n for n in names(df))
    for c in candidates
        if haskey(name_map, c)
            return name_map[c]
        end
    end
    error("Could not find any of these columns: $(candidates). Found: $(names(df))")
end

function parse_coordinate_string(s)
    ss = strip(String(s))
    ss = replace(ss, "(" => ""); ss = replace(ss, ")" => "")
    parts = split(ss, ",")
    length(parts) == 2 || error("Could not parse coordinate string: $s")
    return parse(Float64, strip(parts[1])), parse(Float64, strip(parts[2]))
end

function haversine_km(lat1, lon1, lat2, lon2)
    R  = 6371.0
    φ1 = deg2rad(lat1); λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2); λ2 = deg2rad(lon2)
    dφ = φ2 - φ1; dλ = λ2 - λ1
    a  = sin(dφ/2)^2 + cos(φ1)*cos(φ2)*sin(dλ/2)^2
    return R * 2 * atan(sqrt(a), sqrt(1-a))
end

###############################################################################
# Data loader  (identical to deterministic model)
###############################################################################

function load_data(excel_file::String, parameter_file::String)
    infra = read_sheet(excel_file, "Infrastructure")
    pax   = read_sheet(excel_file, "PassengerGroups")
    # Pool B: on-demand passenger groups, served only in the second stage.
    # Optional sheet; if absent, Pool B is empty and the model is unchanged.
    paxB = try
        read_sheet(excel_file, "PassengerGroupsB")
    catch
        @warn "Sheet 'PassengerGroupsB' not found — on-demand pool B is empty."
        DataFrame()
    end
    plane = read_sheet_any(excel_file, ["PlaneData"])

    id_col   = find_col(infra, [:id])
    pads_col = find_col(infra, [:number_of_parking_pads, :parking_pads, :pads])

    coord_col = any(Symbol(String(n)) == :coordinates for n in names(infra)) ?
                find_col(infra, [:coordinates]) : nothing
    lat_col   = any(Symbol(String(n)) == :latitude  for n in names(infra)) ?
                find_col(infra, [:latitude])  :
                any(Symbol(String(n)) == :lat for n in names(infra)) ?
                find_col(infra, [:lat]) : nothing
    lon_col   = any(Symbol(String(n)) == :longitude for n in names(infra)) ?
                find_col(infra, [:longitude]) :
                any(Symbol(String(n)) == :lon for n in names(infra)) ?
                find_col(infra, [:lon]) : nothing

    group_col = find_col(pax, [:group, :group_id, :id])
    orig_col  = find_col(pax, [:origin])
    dest_col  = find_col(pax, [:destination])
    time_col  = find_col(pax, [:time, :arrival_time])
    q_col     = find_col(pax, [:number_of_passengers, :passengers, :numberofpassengers])
    stop_col  = find_col(pax, [:stopover_allowed, :stopoverallowed, :stopover])

    V = sort(unique(Int.(infra[!, id_col])))

    plane_id_col = find_col(plane, [:plane_id, :id, :evtol_id, :plane])
    base_vp_col  = find_col(plane, [:base_vertiport, :base_vertiport_id, :base, :home_vertiport])
    N = sort(Int.(plane[!, plane_id_col]))
    isempty(N) && error("PlaneData sheet is empty.")
    length(unique(N)) == length(N) || error("Duplicate Plane IDs.")

    bv = Dict{Int,Int}(Int(r[plane_id_col]) => Int(r[base_vp_col]) for r in eachrow(plane))
    bad = sort([b for b in values(bv) if !(b in V)])
    isempty(bad) || error("Invalid Base Vertiport values $(bad). Valid: $(V).")

    A = sort(Int.(pax[!, group_col]))

    # Parameters
    params = Dict{String,Float64}()
    XLSX.openxlsx(parameter_file) do xf
        for row in XLSX.eachrow(xf["Parameters"])
            key, val = row[1], row[2]
            (key !== missing && val !== missing) || continue
            params[String(key)] = parse(Float64, replace(string(val), "," => "."))
        end
    end

    fare_direct_per_km    = params["fare_direct_per_km"]
    fare_stopover_factor  = params["fare_stopover_factor"]
    p_penalty             = params["p_penalty"]
    operating_cost_per_km = params["operating_cost_per_km"]
    battery_per_km        = params["battery_per_km"]
    time_per_km           = params["time_per_km"]
    cap_u                 = params["cap_u"]
    opening_cost          = params["opening_cost"]
    bmax                  = params["bmax"]
    bmid                  = params["bmid"]
    b_penalty             = params["b_penalty"]
    bmin                  = params["bmin"]
    ec                    = params["ec"]
    te                    = params["te"]
    w                     = params["w"]
    ET                    = params["ET"]
    M1  = 1;  M2 = bmax;  M3 = ET
 


    M          = 0:10
    M_no0      = 1:maximum(M)
    M_mid      = 1:(maximum(M)-1)
    M_no_last  = 0:(maximum(M)-1)
    T          = 0:ET
    T_no0      = 1:maximum(T)

    lat   = Dict{Int,Float64}()
    lon   = Dict{Int,Float64}()
    cap_v = Dict{Int,Int}()
    for r in eachrow(infra)
        j = Int(r[id_col])
        cap_v[j] = Int(r[pads_col])
        if coord_col !== nothing
            lat[j], lon[j] = parse_coordinate_string(r[coord_col])
        else
            (lat_col === nothing || lon_col === nothing) &&
                error("Infrastructure needs coordinates or lat/lon columns.")
            lat[j] = Float64(r[lat_col])
            lon[j] = Float64(r[lon_col])
        end
    end

    # Prices
    prices    = read_sheet(parameter_file, "Prices")
    from_col  = find_col(prices, [:from])
    to_col    = find_col(prices, [:to])
    fd_col    = find_col(prices, [:fd_sum])
    dt_col    = find_col(prices, [:bil_tid_min])

    fd_lookup = Dict{Tuple{Int,Int},Float64}()
    drive_time_lookup = Dict{Tuple{Int,Int},Float64}()
    for r in eachrow(prices)
        i, j = Int(r[from_col]), Int(r[to_col])
        fd_lookup[(i,j)]         = Float64(r[fd_col])
        drive_time_lookup[(i,j)] = Float64(r[dt_col])
    end
    for (i,j) in collect(keys(fd_lookup))
        fd_lookup[(j,i)]         = fd_lookup[(i,j)]
        drive_time_lookup[(j,i)] = drive_time_lookup[(i,j)]
    end

    end_vp = Dict{Int,Vector{Int}}()
    for i in V
        end_vp[i] = [j for j in V if j == i ||
                     (haskey(drive_time_lookup,(i,j)) && drive_time_lookup[(i,j)] <= 60.0)]
    end

    dist = Dict{Tuple{Int,Int},Float64}()
    fd   = Dict{Tuple{Int,Int},Float64}()
    fs   = Dict{Tuple{Int,Int},Float64}()
    c    = Dict{Tuple{Int,Int},Float64}()
    e    = Dict{Tuple{Int,Int},Float64}()
    rt   = Dict{Tuple{Int,Int},Int}()
    for i in V, j in V
        dij       = haversine_km(lat[i], lon[i], lat[j], lon[j])
        dist[(i,j)] = dij
        fd[(i,j)]   = get(fd_lookup, (i,j), 0.0)
        fs[(i,j)]   = fare_stopover_factor * fd[(i,j)]
        c[(i,j)]    = operating_cost_per_km * dij
        e[(i,j)]    = battery_per_km * dij
        rt[(i,j)]   = Int(ceil(time_per_km * dij))
    end

    # Passenger parameters
    op = Dict{Int,Int}(); dp = Dict{Int,Int}()
    dt = Dict{Int,Float64}(); q = Dict{Int,Int}()
    so = Dict{Int,Int}(); p = Dict{Int,Float64}()
    for r in eachrow(pax)
        a     = Int(r[group_col])
        op[a] = Int(r[orig_col]);  dp[a] = Int(r[dest_col])
        dt[a] = Float64(r[time_col]); q[a] = Int(r[q_col])
        so[a] = Int(r[stop_col]);  p[a]  = p_penalty
    end

    ###########################################################################
    # Pool B (on-demand) passengers — second stage only.
    # IDs offset by POOL_B_OFFSET so they never collide with Pool A IDs. They
    # share the op/dp/dt/q/so dicts; p[b]=0 (no penalty — opportunistic only).
    ###########################################################################
    POOL_B_OFFSET = 1000
    B = Int[]
    if nrow(paxB) > 0
        gB = find_col(paxB, [:group, :group_id, :id])
        oB = find_col(paxB, [:origin]); dB = find_col(paxB, [:destination])
        tB = find_col(paxB, [:time, :arrival_time])
        qB = find_col(paxB, [:number_of_passengers, :passengers, :numberofpassengers])
        sB = find_col(paxB, [:stopover_allowed, :stopoverallowed, :stopover])
        for r in eachrow(paxB)
            b = POOL_B_OFFSET + Int(r[gB]); push!(B, b)
            op[b] = Int(r[oB]); dp[b] = Int(r[dB]); dt[b] = Float64(r[tB])
            q[b]  = Int(r[qB]); so[b] = Int(r[sB]); p[b]  = 0.0
        end
        B = sort(B)
    end

    # Combined service set: every group that can be SERVED in the second stage.
    AB = vcat(A, B)

    d = Dict((a,i,j) => (i==op[a] && j==dp[a] ? 1 : 0)
             for a in AB, i in V, j in V)

    # Generate scenario-dependent parameters
    rt_s, e_s, S, pi_s = generate_scenarios(V, lat, lon, rt, e, time_per_km)

    ###########################################################################
    # On-demand (Pool B) demand realizations: each Pool B group independently
    # appears in scenario sc with probability rho, fixed by a deterministic seed
    # so the realization matches the heuristic (reproducible, comparable).
    # appear[sc] = set of Pool B group IDs present in scenario sc.
    ###########################################################################
    demand_rho  = 0.5
    demand_seed = 20260101
    appear = Dict{Int, Vector{Int}}()
    for sc in S
        rng = Random.MersenneTwister(demand_seed + sc)
        present = Int[]
        for b in B
            rand(rng) < demand_rho && push!(present, b)
        end
        appear[sc] = sort(present)
    end

    return (
        V=V, A=A, B=B, AB=AB, N=collect(N), S=S, appear=appear,
        M=collect(M), M_no0=collect(M_no0), M_mid=collect(M_mid),
        M_no_last=collect(M_no_last), T=collect(T), T_no0=collect(T_no0),
        bv=bv, lat=lat, lon=lon, dist=dist,
        fd=fd_lookup, fs=fs, c=c,
        e=e,   rt=rt,          # deterministic parameters (kept for reference)
        e_s=e_s, rt_s=rt_s,   # scenario-dependent parameters
        pi_s=pi_s,
        end_vp=end_vp, op=op, dp=dp, dt=dt, q=q, so=so, p=p, d=d,
        cap_v=cap_v, cap_u=cap_u, opening_cost=opening_cost,
        bmax=bmax, bmid=bmid, b_penalty=b_penalty, bmin=bmin,
        ec=ec, te=te, w=w, ET=ET,
        M1=M1, M2=M2, M3=M3,
    )
end

###############################################################################
# Model builder — two-stage stochastic
###############################################################################

function build_model(excel_file::String, parameter_file::String;
                     show_progress::Bool=true, display_interval_sec::Int=5)

    data = load_data(excel_file, parameter_file)

    A=data.A; B=data.B; AB=data.AB; V=data.V; N=data.N; S=data.S
    appear=data.appear
    M=data.M; M_no0=data.M_no0; M_mid=data.M_mid
    M_no_last=data.M_no_last; T=data.T; T_no0=data.T_no0

    bv=data.bv; fd=data.fd; fs=data.fs; c=data.c
    e_s=data.e_s; rt_s=data.rt_s; pi_s=data.pi_s
    d=data.d; op=data.op; dp=data.dp; dt=data.dt
    q=data.q; so=data.so; p=data.p; end_vp=data.end_vp
    cap_v=data.cap_v; cap_u=data.cap_u; opening_cost=data.opening_cost
    bmax=data.bmax; bmid=data.bmid; b_penalty=data.b_penalty; bmin=data.bmin
    ec=data.ec; te=data.te; w=data.w
    M1=data.M1; M2=data.M2; M3=data.M3

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "OutputFlag", show_progress ? 1 : 0)
    show_progress && set_optimizer_attribute(model, "DisplayInterval",
                                             max(1, display_interval_sec))

    ###########################################################################
    # FIRST-STAGE variables  (no scenario index)
    ###########################################################################

    # y[n] = 1 if eVTOL n is activated
    @variable(model, y[n in N], Bin)

    # alpha[a] = 1 if pre-booked passenger group a is ACCEPTED in the first stage.
    # Acceptance is eVTOL-independent: the commitment is to serve the group, not
    # to a specific aircraft. Which eVTOL actually serves an accepted group is a
    # second-stage decision that may differ by scenario.
    @variable(model, alpha[a in A], Bin)

    ###########################################################################
    # SECOND-STAGE variables  (scenario index sc added)
    ###########################################################################

    # x[sc,i,j,m,n] = 1 if eVTOL n flies i->j in operation m under scenario sc
    @variable(model, x[sc in S, i in V, j in V, m in M, n in N], Bin)

    # z[sc,a] = 1 if accepted Pool A group a is routed via a stopover in scenario
    # sc. This is a SECOND-STAGE (operational) decision: whether to route direct
    # or via an intermediate vertiport may differ by scenario, since wind affects
    # which routes are feasible. Only meaningful for Pool A; Pool B is always direct.
    @variable(model, z[sc in S, a in A], Bin)

    # s[sc,a,m,n] = 1 if eVTOL n serves group a in operation m under scenario sc
    # Indexed over AB (Pool A ∪ Pool B): both pools can be served in the second
    # stage. Commitment (ss) and penalties apply to Pool A only; Pool B service
    # is gated by the per-scenario demand realization (fixed to 0 if absent).
    @variable(model, s[sc in S, a in AB, m in M, n in N], Bin)

    # k[sc,a,i,j,m,n] = 1 if eVTOL n flies i->j in op m serving group a, scenario sc
    @variable(model, k[sc in S, a in AB, i in V, j in V, m in M, n in N], Bin)

    # u[sc,m,n] = battery level of eVTOL n after operation m in scenario sc
    @variable(model, u[sc in S, m in M, n in N] >= 0)

    # dep[sc,m,n] / arr[sc,m,n] = departure / arrival time
    @variable(model, dep[sc in S, m in M, n in N] >= 0)
    @variable(model, arr[sc in S, m in M, n in N] >= 0)

    # charge[sc,m,n] = battery charged during turnaround before op m
    @variable(model, charge[sc in S, m in M, n in N] >= 0)

    # over_bmid[sc,m,n] = battery above bmid (penalised in objective).
    # Integer, matching the deterministic model's bp[m,n] variable.
    @variable(model, over_bmid[sc in S, m in M, n in N] >= 0, Int)

    # is_p[sc,j,n,t] = 1 if eVTOL n is parked at j at time t in scenario sc
    @variable(model, is_p[sc in S, j in V, n in N, t in T], Bin)

    # is_o[sc,i,j,m,n,t] = 1 if eVTOL n flies i->j in op m at time t, scenario sc
    @variable(model, is_o[sc in S, i in V, j in V, m in M, n in N, t in T], Bin)

    ###########################################################################
    # Fix operation-0 variables to zero (operation 0 is not a real flight)
    ###########################################################################
    for sc in S, i in V, j in V, n in N
        fix(x[sc,i,j,0,n], 0.0; force=true)
        for a in AB
            fix(k[sc,a,i,j,0,n], 0.0; force=true)
        end
    end
    for sc in S, a in AB, n in N
        fix(s[sc,a,0,n], 0.0; force=true)
    end

    ###########################################################################
    # Demand gating: a Pool B group may be served in scenario sc ONLY if it
    # appears in that scenario's demand realization. If it does not appear, all
    # its service and service-arc variables are fixed to zero in that scenario.
    ###########################################################################
    for sc in S
        present = Set(appear[sc])
        for b in B
            if !(b in present)
                for m in M, n in N
                    fix(s[sc,b,m,n], 0.0; force=true)
                    for i in V, j in V
                        fix(k[sc,b,i,j,m,n], 0.0; force=true)
                    end
                end
            end
        end
    end

    ###########################################################################
    # Objective: maximize expected profit
    #
    # First-stage terms (scenario-independent):
    #   + revenue from accepted passengers
    #   - penalty for unserved passengers
    #   - fleet activation cost
    #
    # Second-stage terms (expected value over scenarios):
    #   - operational routing cost
    #   - battery-above-bmid penalty
    ###########################################################################
    @objective(model, Max,
        # Revenue from accepted prebooked passengers (first-stage), per passenger
        sum(d[(a,i,j)] * alpha[a] * q[a] * (fd[i,j]*(1-so[a]) + fs[i,j]*so[a])
            for a in A, i in V, j in V)
        # Penalty for unaccepted passenger groups (first-stage)
        - sum(p[a] * (1 - alpha[a]) for a in A)
        # Fleet activation cost (first-stage)
        - sum(opening_cost * y[n] for n in N)
        # Expected second-stage value
        + sum(pi_s[sc] * (
            # On-demand (Pool B) revenue: per passenger, earned per scenario when served.
            sum(d[(b,i,j)] * s[sc,b,m,n] * q[b] *
                (fd[i,j]*(1-so[b]) + fs[i,j]*so[b])
                for b in B, i in V, j in V, m in M_no0, n in N)
            # Routing cost under scenario sc
            - sum(c[(i,j)] * x[sc,i,j,m,n]
                  for i in V, j in V, m in M, n in N)
            # Battery penalty under scenario sc
            - sum(over_bmid[sc,m,n] for m in M, n in N) * b_penalty
          )
          for sc in S)
    )

    ###########################################################################
    # NON-ANTICIPATIVITY is enforced structurally: y, ss, z carry no sc index.
    # All sc-indexed second-stage variables must be consistent with these.
    ###########################################################################

    ###########################################################################
    # FIRST-STAGE constraints (no sc index)
    ###########################################################################

    # (6.13) With eVTOL-independent acceptance, the per-eVTOL uniqueness of the
    # old ss[a,n] is no longer needed: constraint (6.7) fixes each accepted
    # group's total number of service legs to alpha[a] + z[sc,a], so a direct
    # accepted group is served on exactly one leg (hence by one eVTOL) and a
    # stopover group on exactly two.
    # (6.11, the stopover-allowed constraint, is now second-stage since z is
    # scenario-dependent — see inside the scenario loop below.)

    ###########################################################################
    # SECOND-STAGE constraints (sc index on all operational variables)
    ###########################################################################

    for sc in S

        # Shorthand accessors for this scenario. The scenario dicts only contain
        # off-diagonal pairs (i != j); a vertiport-to-itself entry has zero travel
        # time and energy, so return 0 on the diagonal rather than a missing key.
        rt_sc = (i,j) -> i == j ? 0   : rt_s[(sc,i,j)]
        e_sc  = (i,j) -> i == j ? 0.0 : e_s[(sc,i,j)]

        # (6.11) Stopover only if passenger allows it (per scenario)
        @constraint(model, [a in A], z[sc,a] <= so[a])

        # ── Fleet activation and routing ─────────────────────────────────────

        # (6.2) eVTOL leaves base vertiport in operation 1 iff activated
        @constraint(model, [n in N],
            sum(x[sc,bv[n],j,1,n] for j in V) == y[n])

        # (6.3a) Return-to-base feasibility for intermediate operations
        @constraint(model, [m in M_mid, n in N],
            sum(x[sc,i,j,m,n] for i in V, j in end_vp[bv[n]]) >=
            sum(x[sc,i,j,m,n] for i in V, j in V) -
            sum(x[sc,i,j,m+1,n] for i in V, j in V))

        # (6.3b) Return-to-base feasibility for last operation
        @constraint(model, [n in N],
            sum(x[sc,i,j,maximum(M),n] for i in V, j in end_vp[bv[n]]) >=
            sum(x[sc,i,j,maximum(M),n] for i in V, j in V))

        # (6.4) At most one route per operation per eVTOL; inactive if y[n]=0
        @constraint(model, [m in M_no0, n in N],
            sum(x[sc,i,j,m,n] for i in V, j in V) <= y[n])

        # (6.5) No self-loops
        @constraint(model, [i in V, m in M_no0, n in N], x[sc,i,i,m,n] <= 0)

        # (6.6) Flow continuity between consecutive operations
        @constraint(model, [j in V, m in M_mid, n in N],
            sum(x[sc,i,j,m,n] for i in V) >= sum(x[sc,j,i2,m+1,n] for i2 in V))

        # ── Passenger assignment ─────────────────────────────────────────────

        # (6.7) Number of legs = acceptance indicator + stopover indicator (Pool A)
        @constraint(model, [a in A],
            sum(s[sc,a,m,n] for m in M_no0, n in N) ==
            alpha[a] + z[sc,a])

        # (6.7-B) Pool B (on-demand): served at most once, directly, and only if
        # it appears in this scenario (appearance enforced by the fix-to-zero
        # gating above). No ss/z — service is purely opportunistic recourse.
        @constraint(model, [b in B],
            sum(s[sc,b,m,n] for m in M_no0, n in N) <= 1)

        # (6.8) Direct connection if z[sc,a]=0 (Pool A)
        @constraint(model, [a in A, m in M, n in N],
            s[sc,a,m,n] <= x[sc,op[a],dp[a],m,n] + z[sc,a]*M1)

        # (6.8-B) Pool B is always direct: served on op[b]->dp[b] leg only.
        @constraint(model, [b in B, m in M, n in N],
            s[sc,b,m,n] <= x[sc,op[b],dp[b],m,n])

        # (6.9) Stopover path upper bound
        @constraint(model, [a in A, m in M, n in N],
            s[sc,a,m,n] <=
            sum(x[sc,op[a],k_node,m,n] for k_node in V) +
            sum(x[sc,k_node,dp[a],m,n] for k_node in V) +
            (1 - z[sc,a])*M1)

        # (6.10a) Stopover path lower bound (m, m+1)
        @constraint(model, [a in A, m in M_mid, n in N],
            s[sc,a,m,n] >=
            sum(x[sc,op[a],k_node,m,n] + x[sc,k_node,dp[a],m+1,n]
                for k_node in V) - 1 - (1-z[sc,a])*M1)

        # (6.10b) Stopover path lower bound (m-1, m)
        @constraint(model, [a in A, m in 2:maximum(M), n in N],
            s[sc,a,m,n] >=
            sum(x[sc,op[a],k_node,m-1,n] + x[sc,k_node,dp[a],m,n]
                for k_node in V) - 1 - (1-z[sc,a])*M1)

        # (6.12) A Pool A group can be served (by any eVTOL, any op) only if it
        # was accepted in the first stage. alpha is eVTOL-independent.
        @constraint(model, [a in A, m in M, n in N],
            s[sc,a,m,n] <= alpha[a])

        # (6.14) At least one arc per served operation (Pool A)
        @constraint(model, [a in A, m in M_no0, n in N],
            sum(k[sc,a,i,j,m,n] for i in V, j in V) >= s[sc,a,m,n])

        # (6.14-B) Same for Pool B
        @constraint(model, [b in B, m in M_no0, n in N],
            sum(k[sc,b,i,j,m,n] for i in V, j in V) >= s[sc,b,m,n])

        # (6.15) Direct service linkage (Pool A)
        @constraint(model, [a in A, i in V, j in V, m in M_no0, n in N],
            2*k[sc,a,i,j,m,n] <= d[(a,i,j)] + x[sc,i,j,m,n] + z[sc,a]*M1)

        # (6.15-B) Pool B direct service linkage (no z; always direct)
        @constraint(model, [b in B, i in V, j in V, m in M_no0, n in N],
            2*k[sc,b,i,j,m,n] <= d[(b,i,j)] + x[sc,i,j,m,n])

        # (6.16) At most 1+z[a] arcs per passenger group (Pool A)
        @constraint(model, [a in A],
            sum(k[sc,a,i,j,m,n] for m in M, i in V, j in V, n in N) <= 1 + z[sc,a])

        # (6.16-B) Pool B uses at most one arc (always direct)
        @constraint(model, [b in B],
            sum(k[sc,b,i,j,m,n] for m in M, i in V, j in V, n in N) <= 1)

        # (6.17a) Stopover service linkage
        @constraint(model, [a in A, i in V, kk in V, j in V, m in M_no_last, n in N],
            k[sc,a,i,kk,m,n] + k[sc,a,kk,j,m+1,n] <=
            x[sc,i,kk,m,n] + x[sc,kk,j,m+1,n] + (1-z[sc,a])*M1)

        # (6.17b)
        @constraint(model, [a in A, i in V, kk in V, j in V, m in M_no_last, n in N],
            k[sc,a,i,kk,m,n] + k[sc,a,kk,j,m+1,n] <= d[(a,i,j)] + 1)

        # (6.18) Seat capacity — sums over BOTH pools (Pool B occupies seats too)
        @constraint(model, [m in M, n in N],
            sum(s[sc,a,m,n]*q[a] for a in AB) <= cap_u)

        # (6.19) Direct-only groups fly alone — over both pools
        @constraint(model, [a in AB, m in M, n in N; so[a]==0],
            sum(s[sc,b,m,n] for b in AB) <= 1 + (length(AB)-1)*(1-s[sc,a,m,n]))

        # ── Battery dynamics ─────────────────────────────────────────────────
        # rt_sc and e_sc now use scenario-dependent parameters,
        # which embed both the directional wind effect (gamma_t) and
        # the temperature effect (phi, via gamma_e).

        # (6.20) eVTOL starts with mid battery
        @constraint(model, [n in N], u[sc,0,n] == bmid)

        # (6.21) Battery cannot exceed max
        @constraint(model, [m in M_no_last, n in N],
            u[sc,m,n] + charge[sc,m+1,n] <= bmax)

        # (6.22) Battery must stay above minimum
        @constraint(model, [m in M, n in N], u[sc,m,n] >= bmin)

        # (6.23) Penalize battery above bmid (includes the energy term, as in
        # the deterministic model). over_bmid plays the role of bp[m,n].
        @constraint(model, [m in M_no0, n in N],
            over_bmid[sc,m,n] >= u[sc,m,n] - bmid +
                sum(e_sc(i,j)*x[sc,i,j,m,n] for i in V, j in V))

        # (6.24a) Battery update between operations
        @constraint(model, [i in V, j in V, m in 1:maximum(M), n in N],
            u[sc,m,n] <= u[sc,m-1,n] - e_sc(i,j)*x[sc,i,j,m,n] +
                         charge[sc,m,n] + (1 - x[sc,i,j,m,n])*M2)

        # (6.24b) Battery update between operations
        @constraint(model, [i in V, j in V, m in 1:maximum(M), n in N],
            u[sc,m,n] >= u[sc,m-1,n] - e_sc(i,j)*x[sc,i,j,m,n] +
                         charge[sc,m,n] - (1 - x[sc,i,j,m,n])*M2)

        # (6.25) Charging level at each operation
        @constraint(model, [i in V, j in V, m in 1:maximum(M), n in N],
            charge[sc,m,n] <= ec*(dep[sc,m,n] - arr[sc,m-1,n]) +
                              (1 - x[sc,i,j,m,n])*M2)

        # ── Timing constraints ───────────────────────────────────────────────
        # All timing constraints use rt_sc (scenario-dependent travel time).

        # (6.26) Operation 0 starts at time 0
        @constraint(model, [n in N], arr[sc,0,n] == 0)

        # (6.27a) Arrival time lower bound for operations >= 2
        # (turnaround time te added because previous flight must have landed
        #  and minimum ground time must elapse before next departure)
        @constraint(model, [m in 2:maximum(M), n in N],
            arr[sc,m,n] >= arr[sc,m-1,n] +
                sum((te + rt_sc(i,j))*x[sc,i,j,m,n] for i in V, j in V))

        # (6.27b) Arrival time for first operation (no prior turnaround)
        @constraint(model, [n in N],
            arr[sc,1,n] >= sum(rt_sc(i,j)*x[sc,i,j,1,n] for i in V, j in V))

        # (6.28) Departure = arrival - travel time
        @constraint(model, [m in M, n in N],
            arr[sc,m,n] == dep[sc,m,n] +
                sum(rt_sc(i,j)*x[sc,i,j,m,n] for i in V, j in V))

        # (6.29) Minimum layover time between consecutive operations (both pools)
        @constraint(model, [a in AB, n in N, m in M_no_last],
            dep[sc,m+1,n] <= arr[sc,m,n] + te +
                             (2 - s[sc,a,m,n] - s[sc,a,m+1,n])*M3)

        # (6.30) Passenger arrives before their eVTOL departs (both pools)
        @constraint(model, [a in AB, i in V, j in V, m in M_no0, n in N],
            d[(a,i,j)]*dt[a] - (1 - (s[sc,a,m,n] - s[sc,a,m-1,n]))*M3 <=
            arr[sc,m,n] - sum(rt_sc(i,kk)*x[sc,i,kk,m,n] for kk in V))

        # (6.31) Maximum passenger waiting time (both pools)
        @constraint(model, [a in AB, m in M_no0, n in N],
            arr[sc,m,n]
            - sum(rt_sc(i,j)*x[sc,i,j,m,n] for i in V, j in V)
            - sum(d[(a,i,j)]*dt[a]*s[sc,a,m,n] for i in V, j in V)
            - (1 - (s[sc,a,m,n] - s[sc,a,m-1,n]))*M3 <= w)

        # ── Time-space occupancy ─────────────────────────────────────────────
        # Occupancy constraints are also scenario-dependent because departure
        # and arrival times differ across scenarios (due to wind-dependent rt).

        # (6.32) eVTOL is either parked or flying at each time t
        @constraint(model, [n in N, t in T],
            sum(is_p[sc,j,n,t] for j in V) +
            sum(is_o[sc,i,j,m,n,t] for i in V, j in V, m in M) <= 1)

        # (6.33) Travel-time occupancy relationship (uses scenario rt).
        # rt_s is already integer-valued, matching the deterministic rt.
        @constraint(model, [i in V, j in V, m in M, n in N],
            rt_sc(i,j) * x[sc,i,j,m,n] ==
            sum(is_o[sc,i,j,m,n,t] for t in T))

        # (6.34) Departure time bound from occupancy
        @constraint(model, [i in V, j in V, m in M_no0, n in N, t in T],
            dep[sc,m,n] <= t + M3*(1 - is_o[sc,i,j,m,n,t]) - 1)

        # (6.35) Arrival time bound from occupancy
        @constraint(model, [i in V, j in V, m in M_no0, n in N, t in T],
            arr[sc,m,n] >= t - M3*(1 - is_o[sc,i,j,m,n,t]))

        # (6.36) Initial parking at base vertiport (tied to activation)
        @constraint(model, [n in N], is_p[sc,bv[n],n,0] == y[n])

        # (6.37) Parking state propagation
        @constraint(model, [j in V, n in N, t in T_no0],
            is_p[sc,j,n,t] <= is_p[sc,j,n,t-1] +
                              sum(is_o[sc,i,j,m,n,t-1] for i in V, m in M))

        # (6.38) Vertiport parking capacity
        @constraint(model, [j in V, t in T],
            sum(is_p[sc,j,n,t] for n in N) <= cap_v[j])

        # (extra, as in deterministic OUR_Model) parking implies some future flight
        @constraint(model, [j in V, n in N, t in T],
            is_p[sc,j,n,t] <=
            sum(is_o[sc,i,k,m,n,t2] for i in V, k in V, m in M, t2 in t:maximum(T)))

    end  # for sc in S

    return model, data
end

###############################################################################
# Solve
###############################################################################

function solve_instance(excel_file::String, parameter_file::String;
                        show_progress::Bool=true, display_interval_sec::Int=5)
    timings = Dict{String,Float64}()

    t_build = @elapsed model, data = build_model(
        excel_file, parameter_file;
        show_progress=show_progress,
        display_interval_sec=display_interval_sec)
    timings["Build model (incl. data load)"] = t_build

    t_opt = @elapsed optimize!(model)
    timings["Optimization"] = t_opt

    term   = termination_status(model)
    primal = primal_status(model)
    println("Termination status: ", term)
    println("Primal status:      ", primal)
    if term in (MOI.OPTIMAL, MOI.TIME_LIMIT, MOI.FEASIBLE_POINT)
        println("Objective value:    ", objective_value(model))
    end

    return model, data, timings
end

function print_timing_summary(timings::Dict{String,Float64})
    println("\n", "="^80, "\nRUNTIME SUMMARY\n", "="^80)
    println(lpad("Program Part", 38), " | ", lpad("Seconds", 10))
    println("-"^53)
    parts = ["Build model (incl. data load)", "Optimization"]
    sub = 0.0
    for part in parts
        haskey(timings, part) || continue
        v = timings[part]; sub += v
        println(lpad(part, 38), " | ", lpad(@sprintf("%.3f", v), 10))
    end
    if haskey(timings, "Total script")
        total = timings["Total script"]
        println("-"^53)
        println(lpad("Measured subtotal", 38), " | ", lpad(@sprintf("%.3f", sub), 10))
        println(lpad("Total script",      38), " | ", lpad(@sprintf("%.3f", total), 10))
    end
end

###############################################################################
# Usage
###############################################################################

excel_file     = joinpath("inputData/inputDataGiant.xlsx")
parameter_file = joinpath("inputData/Parameters.xlsx")
println("Using Excel file: ", excel_file)

total_start = time()
model, data, timings = solve_instance(excel_file, parameter_file)
timings["Total script"] = time() - total_start
print_timing_summary(timings)
