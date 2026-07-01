###############################################################################
# Experiment2_realised.jl
#
# EXPERIMENT 2 — Full-day realised-scenario evaluation.
#
# MULTI-INSTANCE VARIANT: instead of one instance evaluated over several
# seeded runs, this runs ONE instance file per run:
#     run 1 -> INSTANCE_FILES[1]
#     run 2 -> INSTANCE_FILES[2]
#     run 3 -> INSTANCE_FILES[3]
# Each run loads its own instance, solves the FIRST STAGE once to fix the
# commitment, then evaluates REALISED_PER_RUN realised operating days against
# that commitment (each day: draw weather + on-demand demand, solve the second
# stage with the commitment fixed, record profit / objective / service / KPIs).
#
# Season via EVTOL_SEASON (default Sommer):
#     EVTOL_SEASON=Vinter julia Experiment2_realised.jl
#
# Output: experiment2_realised_<season>.csv
###############################################################################

using XLSX, DataFrames, Random, CSV, Printf, Statistics

###############################################################################
# CONFIGURATION
###############################################################################
const SEASON_NAME     = get(ENV, "EVTOL_SEASON", "Sommer")

const SMOKE_TEST = get(ENV, "SMOKE_TEST", "0") == "1"

# One instance file per run. RUNS is derived from the number of files so the two
# can never disagree. Each run r loads INSTANCE_FILES[r].
const INSTANCE_FILES = [
    "LTMexp2Inst1.xlsx",
    "LTMexp2Inst2.xlsx",
    "LTMexp2Inst3.xlsx",
]
const RUNS             = SMOKE_TEST ? length(INSTANCE_FILES) : length(INSTANCE_FILES)
const REALISED_PER_RUN = SMOKE_TEST ? 2 : 10    # realised days evaluated per run
const BASE_FS_SEED    = 12347  # first-stage seed; run r uses BASE_FS_SEED + r
const BASE_REAL_SEED  = 4244   # realised-draw seed; run r uses BASE_REAL_SEED + r

# Two-stage parameters (full quality — this is the headline large-instance run).
const MAXTIME_1ST        = SMOKE_TEST ? Int32(10) : Int32(90)
const N_RESTARTS         = SMOKE_TEST ? 1 : 2
const N_OUTER_ITERS      = SMOKE_TEST ? 2 : 4
const MAXTIME_2ND_SEARCH = SMOKE_TEST ? Int32(3)  : Int32(15)
const MAXTIME_2ND_FINAL  = SMOKE_TEST ? Int32(5)  : Int32(30)
const PRICE_BOOST        = 8.0
const HARD_PENALTY       = 50_000.0
const TOP_C              = 4
const MAXTURN            = 100
const DEMAND_RHO         = 0.5
const DEMAND_SEED        = 20260101

const WIND_VECTORS = Dict(
    "N-Mild" => ( 0.0,  -9.3),  "S-Mild" => ( 0.0,   9.3),
    "E-Mild" => (-9.3,   0.0),  "W-Mild" => ( 9.3,   0.0),
    "N-Mod"  => ( 0.0, -27.0),  "S-Mod"  => ( 0.0,  27.0),
    "E-Mod"  => (-27.0,  0.0),  "W-Mod"  => (27.0,   0.0),
    "N-Str"  => ( 0.0, -46.0),  "S-Str"  => ( 0.0,  46.0),
    "E-Str"  => (-46.0,  0.0),  "W-Str"  => (46.0,   0.0),
)
const PHI_LEVELS = [1.00, 1.20, 1.35]   # AboveZero / Cold / VeryCold

quiet(f) = redirect_stdout(f, devnull)

###############################################################################
# Load sources
###############################################################################
src_dir = joinpath(@__DIR__, "..", "Heuristic", "src")
for f in ["DataLoadFunc.jl","FeasibilityFunc.jl","InitialSolFunc.jl",
          "PassAssignFunc.jl","FitnessFunc.jl","NeighborhoodFunc.jl",
          "SANeighborhood.jl","HeuristicFunc.jl","KPIFunc.jl"]
    include(joinpath(src_dir, f))
end
include(joinpath(@__DIR__, "scenario_generation.jl"))
include(joinpath(@__DIR__, "StochasticHeuristic.jl"))

###############################################################################
# Realised-day sampling
###############################################################################
function sample_index(weights::Vector{Float64}, u::Float64)
    c = 0.0
    for (i, w) in enumerate(weights)
        c += w
        u <= c && return i
    end
    return length(weights)
end

function draw_weather(season::Symbol, rng)
    wind_keys = collect(keys(WIND_SHARES))
    wind_w    = [WIND_SHARES[k] for k in wind_keys]
    w_idx     = sample_index(wind_w ./ sum(wind_w), rand(rng))
    wkey      = wind_keys[w_idx]
    wx, wy    = WIND_VECTORS[wkey]

    tp        = collect(TEMP_PROBS[season])
    t_lvl     = sample_index(tp ./ sum(tp), rand(rng))
    phi       = PHI_LEVELS[t_lvl]
    t_lab     = ("AboveZero", "Cold", "VeryCold")[t_lvl]

    return (wx = wx, wy = wy, phi = phi, label = "$(wkey) / $(t_lab)")
end

function draw_demand(data, rng; rho = 0.5)
    appear = Int[]
    for b in data.B
        rand(rng) < rho && push!(appear, b)
    end
    return sort(appear)
end

function build_realised_day_data(data, weather, B_real, time_per_km)
    v_air = 1.0 / time_per_km
    v_min = 0.1 * v_air
    kmh_to_kmmin = 1.0 / 60.0

    sc_rt   = Dict{Tuple{Int,Int}, Int}()
    sc_e    = Dict{Tuple{Int,Int}, Float64}()
    sc_dist = Dict{Tuple{Int,Int}, Float64}()
    for i in data.V, j in data.V
        if i == j
            sc_rt[(i,i)] = 0; sc_e[(i,i)] = 0.0; sc_dist[(i,i)] = 0.0
        else
            h  = headwind_component(data.lat[i], data.lon[i],
                                    data.lat[j], data.lon[j], weather.wx, weather.wy)
            gt = v_air / max(v_min, v_air + h * kmh_to_kmmin)
            ge = weather.phi * gt
            sc_rt[(i,j)]   = Int(round(gt * data.rt[(i,j)]))
            sc_e[(i,j)]    = ge * data.e[(i,j)]
            sc_dist[(i,j)] = sc_e[(i,j)] / data.battery_per_km
        end
    end

    wind_speed = sqrt(weather.wx^2 + weather.wy^2)
    severe     = (wind_speed >= 45.0) || (weather.phi >= 1.35)
    w_day  = Float64(data.w)  + (severe ? 15.0 : 0.0)
    ET_day = Float64(data.ET) + (severe ? 60.0 : 0.0)

    sc_data = (
        infra=data.infra, pax=data.pax, plane=data.plane,
        V=data.V, A=data.A, B=data.B, N=data.N,
        M=data.M, M_no0=data.M_no0, M_mid=data.M_mid, M_no_last=data.M_no_last,
        T=data.T, T_no0=data.T_no0, bv=data.bv, lat=data.lat, lon=data.lon,
        dist=sc_dist, fd=data.fd, fs=data.fs, c=data.c, e=sc_e, rt=sc_rt,
        end_vp=data.end_vp, op=data.op, dp=data.dp, dt=data.dt, q=data.q,
        so=data.so, p=data.p, d=data.d, cap_v=data.cap_v, cap_u=data.cap_u,
        opening_cost=data.opening_cost, bmax=data.bmax, bmid=data.bmid,
        b_penalty=data.b_penalty, bmin=data.bmin, ec=data.ec, te=data.te,
        w=w_day, ET=ET_day, M1=data.M1, M2a=data.M2a, M2b=data.M2b,
        M2c=data.M2c, M3=data.M3, battery_per_km=data.battery_per_km,
        B_realized=B_real,
    )

    Vmax = maximum(data.V)
    rt_mat = zeros(Int, Vmax, Vmax)
    for (ij,v) in sc_rt
        rt_mat[ij[1], ij[2]] = v
    end
    return sc_data, rt_mat, severe
end

###############################################################################
# Header
###############################################################################
println("="^72)
if SMOKE_TEST
    println("EXPERIMENT 2 — *** SMOKE TEST ***  (one instance per run)")
    println("  This validates the pipeline only; results are NOT for the thesis.")
else
    println("EXPERIMENT 2 — realised full-day evaluation, one instance per run   (season: $SEASON_NAME)")
end
println("="^72)

parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")

out_path = joinpath(@__DIR__,
    SMOKE_TEST ? "experiment2_realised_$(lowercase(SEASON_NAME))_SMOKE.csv"
               : "experiment2_realised_$(lowercase(SEASON_NAME)).csv")
rows = NamedTuple[]

###############################################################################
# One run per instance file.
###############################################################################
for run in 1:RUNS
    inst_file  = INSTANCE_FILES[run]
    excel_file = joinpath(@__DIR__, "..", "inputData", inst_file)
    println("\n--- Run $run/$RUNS : instance $inst_file ---")
    if !isfile(excel_file)
        @warn "Instance file not found, skipping" file=inst_file
        continue
    end

    # Load THIS run's instance and build its planning scenarios.
    data        = load_data(excel_file, parameter_file, excel_file)
    time_per_km = derive_time_per_km(data)
    rt_s, e_s, S, pi_s = generate_scenarios(data.V, data.lat, data.lon,
                                            data.rt, data.e, time_per_km)
    S_active = [sc for sc in collect(S) if pi_s[sc] > 0.0]
    pi_s     = Dict(sc => pi_s[sc] for sc in S_active)
    S        = S_active
    poolB    = Set(data.B)

    # First stage once for this instance, reproducible per run.
    Random.seed!(BASE_FS_SEED + run)
    result = quiet() do
        stochastic_heuristic(
            data, rt_s, e_s, S, pi_s;
            maxTurnaround = MAXTURN, MaxTime_1st = MAXTIME_1ST,
            MaxTime_2nd_search = MAXTIME_2ND_SEARCH,
            MaxTime_2nd_final  = MAXTIME_2ND_FINAL,
            top_c = TOP_C, price_boost = PRICE_BOOST, hard_penalty = HARD_PENALTY,
            n_restarts = N_RESTARTS, n_outer_iters = N_OUTER_ITERS,
            demand_rho = DEMAND_RHO, demand_seed = DEMAND_SEED,
        )
    end
    RP  = result.expected_obj
    fsd = FirstStageDecision(result.accepted_passengers, result.active_evtols,
                             result.tentative_sol, result.first_stage_cost)

    committed = sort(collect(result.accepted_passengers))
    active    = sort(collect(result.active_evtols))
    @printf("  First stage: committed=%s  active=%s  in-sample RP=%+.2f\n",
            string(committed), string(active), RP)

    # Realised days for this run, reproducible per run.
    rng = MersenneTwister(BASE_REAL_SEED + run)
    for rd in 1:REALISED_PER_RUN
        weather = draw_weather(Symbol(SEASON_NAME), rng)
        B_real  = draw_demand(data, rng; rho = DEMAND_RHO)
        sc_data, rt_mat, severe = build_realised_day_data(data, weather, B_real, time_per_km)

        r = quiet() do
            solve_second_stage(fsd, sc_data, rt_mat;
                               maxTurnaround = MAXTURN,
                               MaxTime_2nd   = MAXTIME_2ND_FINAL,
                               top_c = TOP_C, price_boost = PRICE_BOOST,
                               hard_penalty = HARD_PENALTY)
        end

        served = Set(ass.group for ass in r.assignments)
        missA  = count(a -> !(a in served), committed)
        onDemB = count(g -> g in poolB, served)

        k = solution_kpis(r.sol, sc_data, rt_mat)

        # Revenue (true fares x q) and operating cost over flown legs.
        rev = 0.0
        for ass in r.assignments
            a = ass.group; i = sc_data.op[a]; j = sc_data.dp[a]
            rev += sc_data.q[a] * (sc_data.fd[(i,j)] * (1 - sc_data.so[a]) +
                                   sc_data.fs[(i,j)] * sc_data.so[a])
        end
        opcost = 0.0
        for plane in r.sol.planes
            for kk in 1:plane.flightLegs
                opcost += sc_data.c[(Int(plane.route[kk]), Int(plane.route[kk+1]))]
            end
        end

        # Operational profit: income minus operating cost only (no battery
        # penalty, no opening cost, no rejection penalty). Comparable to
        # Experiment 1's actual_profit and used for the heuristic comparison.
        realised_actual = rev - opcost

        # True objective for this realised day, matching the Gurobi @objective.
        # r.profit (from second_stage_profit) already equals the Gurobi
        # second-stage value for this scenario: revenue - routing cost - battery
        # penalty (over_bmid * b_penalty) - miss penalty. Subtracting the
        # first-stage cost (opening cost + rejection penalty, paid once) gives the
        # full objective, exactly as expected_objective assembles it in Exp 1.
        realised_objective = r.profit - result.first_stage_cost

        row = (
            season             = SEASON_NAME,
            run                = run,
            instance           = inst_file,
            realised_idx       = rd,
            weather_label      = weather.label,
            severe             = severe ? 1 : 0,
            in_sample_RP       = RP,
            realised_objective = realised_objective,
            realised_actual_profit = realised_actual,
            solver_profit      = r.profit,
            n_committed_poolA  = length(committed),
            n_active_evtols    = length(active),
            first_stage_cost   = result.first_stage_cost,
            committed_missed   = missA,
            poolB_appeared     = length(B_real),
            poolB_served       = onDemB,
            aircraft_utilization           = k.aircraft_utilization,
            deadhead_time                  = k.deadhead_time,
            passenger_demand_served        = k.passenger_demand_served,
            average_passenger_waiting_time = k.average_passenger_waiting_time,
        )
        push!(rows, row)

        @printf("    day %d/%d : %-22s  obj=%+.2f  actual=%+.2f  missA=%d  Bserved=%d/%d  %s\n",
                rd, REALISED_PER_RUN, weather.label, realised_objective, realised_actual, missA,
                onDemB, length(B_real), severe ? "(severe)" : "")

        CSV.write(out_path, DataFrame(rows))   # incremental save
    end
end

println("\n" * "="^72)
println("Done. $(length(rows)) realised-day rows written to:")
println("  $out_path")
println("="^72)