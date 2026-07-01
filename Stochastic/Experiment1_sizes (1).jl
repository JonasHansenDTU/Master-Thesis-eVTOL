###############################################################################
# Experiment1_sizes.jl
#
# EXPERIMENT 1 — Two-stage stochastic heuristic across the six instance sizes
# (Table 9.1), so its results can be placed alongside the existing exact-vs-
# deterministic comparison (Table 9.2).
#
# For each instance size, the two-stage heuristic is run several times (reps).
# Each rep records:
#   - the stochastic objective E[total profit] (RP)            -> "H.Obj"
#   - the wall-clock runtime                                   -> "H.Time"
#   - the first-stage commitment (n committed Pool A, n active eVTOLs)
#   - probability-weighted KPIs across the planning scenarios:
#       aircraft_utilization, deadhead_time,
#       passenger_demand_served, average_passenger_waiting_time
#
# All rows (every rep, not just averages) are appended to a CSV so the thesis
# can report mean and standard deviation per size. Compute the summary in the
# write-up (or a small post-processing step) from the raw rows.
#
# SEASON is selected via the EVTOL_SEASON environment variable so summer and
# winter run as separate processes without editing any source file:
#
#     EVTOL_SEASON=Sommer julia Experiment1_sizes.jl
#     EVTOL_SEASON=Vinter julia Experiment1_sizes.jl
#
# Output: experiment1_<season>.csv  (one file per season)
###############################################################################

using XLSX, DataFrames, Random, CSV, Printf, Statistics

###############################################################################
# CONFIGURATION
###############################################################################

# Season comes from the environment (default Sommer). scenario_generation.jl
# reads the same variable, so the planning scenarios match.
const SEASON_NAME = get(ENV, "EVTOL_SEASON", "Sommer")

# ── SMOKE TEST ───────────────────────────────────────────────────────────────
# Set the environment variable SMOKE_TEST=1 to run a tiny, fast validation:
# 1 rep on only the two smallest instances, with minimal time budgets. This
# confirms the full pipeline (load → solve → KPIs → CSV columns) works in a
# couple of minutes BEFORE committing to the multi-hour real run.
#     SMOKE_TEST=1 julia Experiment1_sizes.jl
# Output goes to experiment1_<season>_SMOKE.csv so it never overwrites real runs.
const SMOKE_TEST = get(ENV, "SMOKE_TEST", "0") == "1"

# Reps per instance. One rep: the heuristic was observed to converge to the same
# solution across reps (identical RP and routing), so a single rep captures the
# result and halves the runtime relative to two reps.
reps_for(id::Int) = 1

# A fixed base seed; each rep uses base + rep so reps differ but the whole
# experiment is reproducible. (Set to a fixed value for a reproducible table.)
const BASE_SEED = 20260601

# The six instances (Table 9.1). Order = scenario ID 1..6.
# (scenario_id, filename, n_VP, n_eVTOL, n_passengers) — the last three are the
# nominal Table-9.1 sizes, recorded in the CSV for the reader's reference.
const ALL_INSTANCES = [
    (id=1, file="inputDataEx2_1_2.xlsx",     vp=2,  evtol=1,  pax=2),
    (id=2, file="inputDataEx5_3_10.xlsx",    vp=5,  evtol=3,  pax=10),
    (id=3, file="inputDataEx5_5_20.xlsx",    vp=5,  evtol=5,  pax=20),
    (id=4, file="inputDataEx5_10_25.xlsx",   vp=5,  evtol=10, pax=25),
    (id=5, file="inputDataEx10_15_35.xlsx",  vp=10, evtol=15, pax=35),
    (id=6, file="inputDataEx10_20_40.xlsx",  vp=10, evtol=20, pax=40),
]

# Smoke test uses only the two smallest instances; the real run uses all six.
const INSTANCES = SMOKE_TEST ? ALL_INSTANCES[1:2] : ALL_INSTANCES

# Two-stage parameters. With several days available, these are set for good
# solution quality. Budgets scale with instance size; the largest instances get
# somewhat reduced per-scenario budgets to keep total runtime within the window,
# but enough to find good commitments. Documented as a size-dependent budget.
function params_for(id::Int)
    if SMOKE_TEST
        return (MaxTime_1st=Int32(10), n_restarts=1, n_outer_iters=2,
                MaxTime_2nd_search=Int32(3), MaxTime_2nd_final=Int32(5))
    elseif id <= 3
        return (MaxTime_1st=Int32(90),  n_restarts=2, n_outer_iters=8,
                MaxTime_2nd_search=Int32(15), MaxTime_2nd_final=Int32(30))
    elseif id <= 4
        return (MaxTime_1st=Int32(90),  n_restarts=2, n_outer_iters=6,
                MaxTime_2nd_search=Int32(12), MaxTime_2nd_final=Int32(25))
    else
        return (MaxTime_1st=Int32(60),  n_restarts=2, n_outer_iters=5,
                MaxTime_2nd_search=Int32(10), MaxTime_2nd_final=Int32(20))
    end
end

# Optionally skip instances whose results you already have (set via START_FROM).
const START_FROM = parse(Int, get(ENV, "START_FROM", "1"))

# Optionally run ONLY a single instance by id (set via ONLY). 0 = disabled (run
# all / from START_FROM). e.g. ONLY=3 runs just instance 3.
const ONLY = parse(Int, get(ENV, "ONLY", "0"))

const PRICE_BOOST  = 8.0
const HARD_PENALTY = 50_000.0
const TOP_C        = 4
const MAXTURN      = 100
const DEMAND_RHO   = 0.5
const DEMAND_SEED  = 20260101

###############################################################################
# Quiet wrapper (suppress the heuristic's verbose internal logging).
###############################################################################
quiet(f) = redirect_stdout(f, devnull)

###############################################################################
# Load sources. (Order matters: data/feasibility/etc. before the heuristic.)
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
# Helper: rebuild a per-scenario integer rt matrix from a scenario's rt dict.
# (scenario_results[sc].sc_data.rt is a Dict; solution_kpis wants a matrix.)
###############################################################################
function rt_matrix_from(sc_data)
    Vmax = maximum(sc_data.V)
    mat  = zeros(Int, Vmax, Vmax)
    for (ij, v) in sc_data.rt
        mat[ij[1], ij[2]] = Int(round(v))
    end
    return mat
end

###############################################################################
# Helper: probability-weighted KPIs across the planning scenarios for one
# heuristic result. Returns a NamedTuple of the four KPIs (weighted means).
###############################################################################
function weighted_kpis(result)
    sr  = result.scenario_results
    pis = result.pi_s
    # Normalise weights over the scenarios actually present in the results.
    scs = collect(keys(sr))
    wsum = sum(pis[sc] for sc in scs)

    util = 0.0; dead = 0.0; served = 0.0; wait = 0.0
    for sc in scs
        r      = sr[sc]
        w      = pis[sc] / wsum
        rt_mat = rt_matrix_from(r.sc_data)
        k      = solution_kpis(r.sol, r.sc_data, rt_mat)
        util   += w * k.aircraft_utilization
        dead   += w * k.deadhead_time
        served += w * k.passenger_demand_served
        wait   += w * k.average_passenger_waiting_time
    end
    return (aircraft_utilization = util,
            deadhead_time = dead,
            passenger_demand_served = served,
            average_passenger_waiting_time = wait)
end

###############################################################################
# Helper: probability-weighted ACTUAL profit (no penalties).
#
# The objective RP (result.expected_obj) includes the commitment hard-penalty
# and the physical feasibility penalties. The "actual profit" reported here is
# the same quantity with those penalties stripped: it is the true economic
# profit the operator would realise, i.e.
#     actual_RP = -first_stage_cost
#                 + sum_sc pi_s * ( revenue_sc - operating_cost_sc )
# computed with the SAME revenue/cost conventions as second_stage_profit
# (true fares x group size, operating cost over flown legs). Opening cost is a
# first-stage cost counted once, not per scenario.
###############################################################################
function weighted_actual_profit(result)
    sr  = result.scenario_results
    pis = result.pi_s
    scs = collect(keys(sr))
    wsum = sum(pis[sc] for sc in scs)

    exp_actual_2nd = 0.0
    for sc in scs
        r  = sr[sc]
        sd = r.sc_data
        w  = pis[sc] / wsum

        # Revenue from served passengers (true fares x q), matching
        # second_stage_profit. r.assignments lists the served groups.
        rev = 0.0
        for ass in r.assignments
            a = ass.group
            i = sd.op[a]; j = sd.dp[a]
            rev += sd.q[a] * (sd.fd[(i,j)] * (1 - sd.so[a]) +
                              sd.fs[(i,j)] * sd.so[a])
        end

        # Operating cost over all flown legs.
        opcost = 0.0
        for plane in r.sol.planes
            for k in 1:plane.flightLegs
                opcost += sd.c[(Int(plane.route[k]), Int(plane.route[k+1]))]
            end
        end

        exp_actual_2nd += w * (rev - opcost)
    end

    # First-stage opening cost is paid once (not per scenario).
    return exp_actual_2nd - result.first_stage_cost
end

###############################################################################
# Run the experiment.
###############################################################################
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")
out_path       = joinpath(@__DIR__,
    SMOKE_TEST ? "experiment1_$(lowercase(SEASON_NAME))_SMOKE.csv"
    : ONLY != 0 ? "experiment1_$(lowercase(SEASON_NAME))_inst$(ONLY).csv"
               : "experiment1_$(lowercase(SEASON_NAME)).csv")

rows = NamedTuple[]

println("="^72)
if SMOKE_TEST
    println("EXPERIMENT 1 — *** SMOKE TEST *** (tiny budgets, 2 instances, 1 rep)")
    println("  This validates the pipeline only; results are NOT for the thesis.")
else
    println("EXPERIMENT 1 — two-stage heuristic across sizes   (season: $SEASON_NAME)")
end
println("="^72)

for inst in INSTANCES
    ONLY != 0 && inst.id != ONLY && continue   # ONLY set: run just that instance
    inst.id < START_FROM && continue   # skip instances already completed
    excel_file = joinpath(@__DIR__, "..", "inputData", "Experiments", inst.file)
    if !isfile(excel_file)
        @warn "Instance file not found, skipping" file=inst.file
        continue
    end

    println("\n--- Instance $(inst.id): $(inst.file)  (VP=$(inst.vp), eVTOL=$(inst.evtol), pax=$(inst.pax)) ---")

    # Load data and build planning scenarios once per instance.
    # Local load_data signature is (infra_file, parameter_file, excel_file).
    # Each instance file is self-contained (Infrastructure + PlaneData +
    # PassengerGroups + PassengerGroupsB), so it serves as both infra and data.
    data        = load_data(excel_file, parameter_file, excel_file)
    time_per_km = derive_time_per_km(data)
    rt_s, e_s, S, pi_s = generate_scenarios(data.V, data.lat, data.lon,
                                            data.rt, data.e, time_per_km)
    # Keep only scenarios with positive probability (zeroed-out by season).
    S_active = [sc for sc in collect(S) if pi_s[sc] > 0.0]
    pi_s     = Dict(sc => pi_s[sc] for sc in S_active)
    S        = S_active

    p = params_for(inst.id)
    n_reps = reps_for(inst.id)

    for rep in 1:n_reps
        Random.seed!(BASE_SEED + rep)   # reproducible, differs per rep

        t0 = time()
        result = quiet() do
            stochastic_heuristic(
                data, rt_s, e_s, S, pi_s;
                maxTurnaround = MAXTURN, MaxTime_1st = p.MaxTime_1st,
                MaxTime_2nd_search = p.MaxTime_2nd_search,
                MaxTime_2nd_final  = p.MaxTime_2nd_final,
                top_c = TOP_C, price_boost = PRICE_BOOST,
                hard_penalty = HARD_PENALTY,
                n_restarts = p.n_restarts, n_outer_iters = p.n_outer_iters,
                demand_rho = DEMAND_RHO, demand_seed = DEMAND_SEED,
            )
        end
        runtime = time() - t0

        k = weighted_kpis(result)
        actual_RP = weighted_actual_profit(result)
        total_unserved = sum(r.n_unserved for r in values(result.scenario_results))

        # Probability-weighted mean number of groups served per scenario, split
        # by pool. Pool B (on-demand) groups carry the +1000 id offset.
        let sr = result.scenario_results, pis = result.pi_s
            scs = collect(keys(sr)); wsum = sum(pis[sc] for sc in scs)
            meanA = 0.0; meanB = 0.0
            for sc in scs
                served = Set(ass.group for ass in sr[sc].assignments)
                nA = count(g -> g < 1000, served)
                nB = count(g -> g >= 1000, served)
                meanA += (pis[sc]/wsum) * nA
                meanB += (pis[sc]/wsum) * nB
            end
            global _mean_poolA_served = meanA
            global _mean_poolB_served = meanB
        end

        row = (
            season              = SEASON_NAME,
            scenario_id         = inst.id,
            n_VP                = inst.vp,
            n_eVTOL             = inst.evtol,
            n_passengers        = inst.pax,
            rep                 = rep,
            H_obj_RP            = result.expected_obj,
            actual_profit_RP    = actual_RP,
            H_time_sec          = runtime,
            n_committed_poolA   = length(result.accepted_passengers),
            n_active_evtols     = length(result.active_evtols),
            first_stage_cost    = result.first_stage_cost,
            committed_misses    = total_unserved,
            mean_poolA_served   = _mean_poolA_served,
            mean_poolB_served   = _mean_poolB_served,
            aircraft_utilization           = k.aircraft_utilization,
            deadhead_time                  = k.deadhead_time,
            passenger_demand_served        = k.passenger_demand_served,
            average_passenger_waiting_time = k.average_passenger_waiting_time,
        )
        push!(rows, row)

        @printf("  rep %d/%d : RP=%+.2f  actual=%+.2f  time=%.1fs  committed=%d  active=%d  misses=%d  poolB_served=%.2f\n",
                rep, n_reps, row.H_obj_RP, row.actual_profit_RP, row.H_time_sec,
                row.n_committed_poolA, row.n_active_evtols, row.committed_misses,
                row.mean_poolB_served)

        # Write incrementally so a crash mid-experiment doesn't lose everything.
        CSV.write(out_path, DataFrame(rows))
    end
end

println("\n" * "="^72)
println("Done. $(length(rows)) rows written to:")
println("  $out_path")
println("="^72)
