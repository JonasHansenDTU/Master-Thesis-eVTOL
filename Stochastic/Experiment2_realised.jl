###############################################################################
# Experiment2_realised.jl
#
# EXPERIMENT 2 — Full-day realised-scenario evaluation on the large instance.
#
# Structure (as discussed):
#   * RUNS independent runs (default 5). Each run:
#       1. Seeds the RNG (per-run) and runs the FIRST STAGE once with n_restarts
#          (exactly as MainCall does) to fix (committed Pool A, active eVTOLs).
#       2. Draws REALISED_PER_RUN realised operating days (default 5). Each day:
#            - draw a realised weather (wind + temperature) from the season dist;
#            - draw a realised on-demand (Pool B) set at rho;
#            - build that day's weather-adjusted data (+ severe slack);
#            - solve the SECOND STAGE once with the commitment fixed (this is
#              ONE realised scenario, not all planning scenarios);
#            - record realised profit, service, and KPIs.
#   * Every realised day (RUNS * REALISED_PER_RUN rows) is saved to a CSV.
#
# This reuses the realised-day machinery (draw_weather / draw_demand /
# build_realised_day_data) and the current solve_second_stage, so it inherits
# the active-fleet confinement, Pool B handling, and ride-sharing fixes.
#
# Season via EVTOL_SEASON (default Sommer):
#     EVTOL_SEASON=Vinter julia Experiment2_realised.jl
#
# NOTE: This uses the large instance (inputDataGiant.xlsx). If you generate a
# new scenario_generation.jl per run, re-run this script after swapping it in;
# alternatively the per-run RNG seed already varies the realised draws each run.
#
# Output: experiment2_realised_<season>.csv
###############################################################################

using XLSX, DataFrames, Random, CSV, Printf, Statistics

###############################################################################
# CONFIGURATION
###############################################################################
const SEASON_NAME     = get(ENV, "EVTOL_SEASON", "Sommer")

# ── SMOKE TEST ───────────────────────────────────────────────────────────────
# SMOKE_TEST=1 runs a tiny, fast validation: 2 runs x 2 realised days with
# minimal budgets, confirming the pipeline and CSV columns before the real run.
#     SMOKE_TEST=1 julia Experiment2_realised.jl
const SMOKE_TEST = get(ENV, "SMOKE_TEST", "0") == "1"

const RUNS            = SMOKE_TEST ? 2 : 5      # independent first-stage runs
const REALISED_PER_RUN = SMOKE_TEST ? 2 : 5     # realised days evaluated per run
const BASE_FS_SEED    = 12345  # first-stage seed; run r uses BASE_FS_SEED + r
const BASE_REAL_SEED  = 4242   # realised-draw seed; run r uses BASE_REAL_SEED + r

# Large instance.
const INSTANCE_FILE = "inputDataGiant.xlsx"

# Two-stage parameters (full quality — this is the headline large-instance run).
# The smoke test overrides these with tiny budgets just below.
const MAXTIME_1ST        = SMOKE_TEST ? Int32(10) : Int32(90)
const N_RESTARTS         = SMOKE_TEST ? 1 : 2
const N_OUTER_ITERS      = SMOKE_TEST ? 2 : 8
const MAXTIME_2ND_SEARCH = SMOKE_TEST ? Int32(3)  : Int32(15)
const MAXTIME_2ND_FINAL  = SMOKE_TEST ? Int32(5)  : Int32(30)
const PRICE_BOOST        = 8.0
const HARD_PENALTY       = 50_000.0
const TOP_C              = 4
const MAXTURN            = 100
const DEMAND_RHO         = 0.5
const DEMAND_SEED        = 20260101

# Temperature levels / seasonal mix (DMI-DATA.xlsx Sheet2) and wind options
# (DMI TR12-19). Identical to Realised_scenario.jl.
const PHI_MAP = [1.00, 1.20, 1.35]
const TEMP_PROBS = Dict(
    :Foraar   => (0.825, 0.153, 0.022),
    :Sommer   => (1.000, 0.000, 0.000),
    :Efteraar => (0.934, 0.055, 0.011),
    :Vinter   => (0.629, 0.293, 0.078),
)
const WIND_OPTIONS = [
    (wx=  0.0, wy=  0.0, key="Calm",        share=0.05),
    (wx=-21.0, wy= 21.0, key="SW_moderate", share=0.15),
    (wx=-39.0, wy= 39.0, key="SW_strong",   share=0.08),
    (wx=-30.0, wy=  0.0, key="W_moderate",  share=0.14),
    (wx=-55.0, wy=  0.0, key="W_strong",    share=0.06),
    (wx= 30.0, wy=  0.0, key="E_moderate",  share=0.15),
    (wx=  0.0, wy= 30.0, key="S_moderate",  share=0.20),
    (wx=  0.0, wy=-30.0, key="N_moderate",  share=0.08),
    (wx=  0.0, wy=-30.0, key="N_strong",    share=0.02),
]

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
# Realised-day sampling (mirrors Realised_scenario.jl)
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
    tp    = collect(TEMP_PROBS[season])
    t_lvl = sample_index(tp ./ sum(tp), rand(rng))
    phi   = PHI_MAP[t_lvl]
    t_lab = ["AboveZero","Cold","VeryCold"][t_lvl]
    wsh   = [w.share for w in WIND_OPTIONS]
    w_idx = sample_index(wsh ./ sum(wsh), rand(rng))
    w     = WIND_OPTIONS[w_idx]
    return (wx=w.wx, wy=w.wy, phi=phi, label="$(w.key) / $(t_lab)")
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
    w_day  = Float64(data.w)  + (severe ? 10.0 : 0.0)
    ET_day = Float64(data.ET) + (severe ? 30.0 : 0.0)

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
# Load data and planning scenarios once.
###############################################################################
println("="^72)
if SMOKE_TEST
    println("EXPERIMENT 2 — *** SMOKE TEST *** (tiny budgets, 2 runs x 2 days)")
    println("  This validates the pipeline only; results are NOT for the thesis.")
else
    println("EXPERIMENT 2 — realised full-day evaluation   (season: $SEASON_NAME)")
end
println("="^72)

excel_file     = joinpath(@__DIR__, "..", "inputData", INSTANCE_FILE)
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")
data        = load_data(excel_file, parameter_file)
time_per_km = derive_time_per_km(data)

rt_s, e_s, S, pi_s = generate_scenarios(data.V, data.lat, data.lon,
                                        data.rt, data.e, time_per_km)
S_active = [sc for sc in collect(S) if pi_s[sc] > 0.0]
pi_s     = Dict(sc => pi_s[sc] for sc in S_active)
S        = S_active

out_path = joinpath(@__DIR__,
    SMOKE_TEST ? "experiment2_realised_$(lowercase(SEASON_NAME))_SMOKE.csv"
               : "experiment2_realised_$(lowercase(SEASON_NAME)).csv")
rows = NamedTuple[]
poolB = Set(data.B)

###############################################################################
# RUNS independent runs; each fixes a commitment then evaluates realised days.
###############################################################################
for run in 1:RUNS
    println("\n--- Run $run/$RUNS ---")

    # First stage once (with n_restarts), reproducible per run.
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

        # KPIs for this realised day's second-stage solution.
        k = solution_kpis(r.sol, sc_data, rt_mat)

        row = (
            season             = SEASON_NAME,
            run                = run,
            realised_idx       = rd,
            weather_label      = weather.label,
            severe             = severe ? 1 : 0,
            in_sample_RP       = RP,
            realised_profit    = r.profit,
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

        @printf("    day %d/%d : %-22s  profit=%+.2f  missA=%d  Bserved=%d/%d  %s\n",
                rd, REALISED_PER_RUN, weather.label, r.profit, missA,
                onDemB, length(B_real), severe ? "(severe)" : "")

        CSV.write(out_path, DataFrame(rows))   # incremental save
    end
end

println("\n" * "="^72)
println("Done. $(length(rows)) realised-day rows written to:")
println("  $out_path")
println("="^72)
