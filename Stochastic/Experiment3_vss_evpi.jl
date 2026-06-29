###############################################################################
# Experiment3_vss_evpi.jl
#
# Heuristic-based estimates of the Value of the Stochastic Solution (VSS) and
# the Expected Value of Perfect Information (EVPI), following Birge & Louveaux,
# "Introduction to Stochastic Programming", Chapter 4 — adapted to a MAXIMISATION
# (profit) problem, so the inequalities are reversed relative to the book's
# minimisation form:
#
#   WS  (wait-and-see)        : E_s[ max over a free first stage in scenario s ]
#   RP  (recourse problem)    : the two-stage stochastic solution (here-and-now)
#   EEV (expected result of   : take the EV (mean-scenario) first-stage decision,
#        the EV solution)        fix it, evaluate it across the real scenarios
#
#   EVPI = WS  - RP   (>= 0)  : value of perfect information
#   VSS  = RP  - EEV  (>= 0)  : value of using the stochastic solution instead
#                               of the naive mean-value (expected-value) solution
#
# Ordering (maximisation):  WS >= RP >= EEV.  Both metrics should be >= 0; a
# negative value signals a computation or search-quality issue.
#
# All three quantities are computed with the SAME solver (the heuristic), so the
# comparison is internally consistent; the results are heuristic ESTIMATES of
# VSS/EVPI, not the exact MIP values, and should be reported as such.
#
# The expected-value (mean) scenario uses probability-weighted average travel
# times and energy, and commits on the pre-booked (Pool A) demand only — Pool B
# is opportunistic on-demand traffic that a mean-weather planner would not
# commit a fleet around. Pool B still appears normally in the per-scenario
# evaluation (EEV) and in RP/WS.
#
# Season via EVTOL_SEASON (default Sommer):
#     EVTOL_SEASON=Sommer julia Experiment3_vss_evpi.jl
#
# Output: vss_evpi_<season>.csv  (one row per instance)
###############################################################################

using XLSX, DataFrames, Random, CSV, Printf, Statistics

const SEASON_NAME = get(ENV, "EVTOL_SEASON", "Sommer")
const SMOKE_TEST  = get(ENV, "SMOKE_TEST", "0") == "1"
const BASE_SEED   = 20260601

# Instances (Table 9.1). VSS/EVPI are reported across sizes to show the trend.
const ALL_INSTANCES = [
    (id=1, file="inputDataEx2_1_2.xlsx",     vp=2,  evtol=1,  pax=2),
    (id=2, file="inputDataEx5_3_10.xlsx",    vp=5,  evtol=3,  pax=10),
    (id=3, file="inputDataEx5_5_20.xlsx",    vp=5,  evtol=5,  pax=20),
    (id=4, file="inputDataEx5_10_25.xlsx",   vp=5,  evtol=10, pax=25),
    (id=5, file="inputDataEx10_15_35.xlsx",  vp=10, evtol=15, pax=35),
    (id=6, file="inputDataEx10_20_40.xlsx",  vp=10, evtol=20, pax=40),
]
# VSS/EVPI is reported on instances 1–4. The two largest instances are excluded
# here because the wait-and-see step solves one heuristic run per scenario, which
# is prohibitively expensive on the large instances; 1–4 give a complete and
# standard VSS/EVPI result. (Run 5–6 separately later if time allows.)
#
# START_FROM lets you resume after a crash: START_FROM=2 skips instances with
# id < 2 (already computed and saved) and begins at instance 2. New rows are
# APPENDED to the existing CSV so earlier instances' results are preserved.
#     START_FROM=2 EVTOL_SEASON=Sommer julia Experiment3_vss_evpi.jl
const START_FROM = parse(Int, get(ENV, "START_FROM", "1"))
const INSTANCES = filter(x -> x.id >= START_FROM,
                         SMOKE_TEST ? ALL_INSTANCES[1:2] : ALL_INSTANCES[1:4])

# Tiered budgets, matching Experiment 1 so the estimates are comparable.
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

const PRICE_BOOST  = 8.0
const HARD_PENALTY = 50_000.0
const TOP_C        = 4
const MAXTURN      = 100
const DEMAND_RHO   = 0.5
const DEMAND_SEED  = 20260101

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
# Helpers
###############################################################################

# Penalty-free probability-weighted actual profit of a heuristic result, using
# TRUE fares (matching Experiment 1's actual_profit). This is the quantity used
# for RP and EEV so all three numbers are on the same (true-profit) footing.
function actual_profit_of(result)
    sr  = result.scenario_results
    pis = result.pi_s
    scs = collect(keys(sr))
    wsum = sum(pis[sc] for sc in scs)
    exp2nd = 0.0
    for sc in scs
        r = sr[sc]; sd = r.sc_data; w = pis[sc]/wsum
        rev = 0.0
        for ass in r.assignments
            a = ass.group; i = sd.op[a]; j = sd.dp[a]
            rev += sd.q[a]*(sd.fd[(i,j)]*(1-sd.so[a]) + sd.fs[(i,j)]*sd.so[a])
        end
        opcost = 0.0
        for plane in r.sol.planes
            for k in 1:plane.flightLegs
                opcost += sd.c[(Int(plane.route[k]), Int(plane.route[k+1]))]
            end
        end
        exp2nd += w*(rev - opcost)
    end
    return exp2nd - result.first_stage_cost
end

# Build a single mean (expected-value) scenario: probability-weighted average of
# rt_s and e_s over the real scenarios. rt_s and e_s are flat dicts keyed by the
# 3-tuple (sc, i, j). The mean scenario reuses scenario id 1 (a valid index into
# the global SCENARIOS weather table, which the heuristic uses only for a printed
# label); the actual travel-time and energy data passed here are the mean values,
# so the weather label is cosmetic and irrelevant.
function build_mean_scenario(rt_s, e_s, S, pi_s)
    wsum = sum(pi_s[sc] for sc in S)
    s0 = first(S)
    ij_pairs = [(i, j) for (sc, i, j) in keys(rt_s) if sc == s0]

    mid = 1   # reuse a valid scenario id; data below overrides the weather
    rt_mean = Dict{Tuple{Int,Int,Int},Float64}()
    e_mean  = Dict{Tuple{Int,Int,Int},Float64}()
    for (i, j) in ij_pairs
        rt_mean[(mid, i, j)] = sum(pi_s[sc] * rt_s[(sc, i, j)] for sc in S) / wsum
        e_mean[(mid, i, j)]  = sum(pi_s[sc] * e_s[(sc, i, j)]  for sc in S) / wsum
    end
    return rt_mean, e_mean, [mid], Dict(mid => 1.0)
end

###############################################################################
# Run
###############################################################################
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")
out_path = joinpath(@__DIR__,
    SMOKE_TEST ? "vss_evpi_$(lowercase(SEASON_NAME))_SMOKE.csv"
               : "vss_evpi_$(lowercase(SEASON_NAME)).csv")
rows = NamedTuple[]

println("="^72)
println("VSS / EVPI (heuristic-based)   season: $SEASON_NAME",
        SMOKE_TEST ? "   *** SMOKE TEST ***" : "")
println("="^72)

for inst in INSTANCES
    excel_file = joinpath(@__DIR__, "..", "inputData", "Experiments", inst.file)
    isfile(excel_file) || (@warn "missing instance, skipping" file=inst.file; continue)

    println("\n--- Instance $(inst.id): $(inst.file) ---")
    data        = load_data(excel_file, parameter_file, excel_file)
    time_per_km = derive_time_per_km(data)
    rt_s, e_s, S, pi_s = generate_scenarios(data.V, data.lat, data.lon,
                                            data.rt, data.e, time_per_km)
    S = [sc for sc in collect(S) if pi_s[sc] > 0.0]
    pi_s = Dict(sc => pi_s[sc] for sc in S)

    p = params_for(inst.id)
    # Single-scenario budget for the WS and EV solves (each is an easy one-
    # scenario problem).
    common = (maxTurnaround=MAXTURN, MaxTime_1st=p.MaxTime_1st,
              MaxTime_2nd_search=p.MaxTime_2nd_search,
              MaxTime_2nd_final=p.MaxTime_2nd_final, top_c=TOP_C,
              price_boost=PRICE_BOOST, hard_penalty=HARD_PENALTY,
              n_restarts=p.n_restarts, n_outer_iters=p.n_outer_iters,
              demand_rho=DEMAND_RHO, demand_seed=DEMAND_SEED)
    # The RP problem spans all scenarios and is much harder; give it extra outer
    # iterations so it is solved to a quality comparable to the easy WS/EV solves.
    # Without this, RP is systematically under-solved and the required ordering
    # WS >= RP >= EEV (for maximisation) can be violated by heuristic noise.
    rp_common = merge(common, (n_outer_iters = p.n_outer_iters + 4,))

    # ----- RP : the stochastic (here-and-now) solution -----------------------
    Random.seed!(BASE_SEED)
    rp_res = quiet() do
        stochastic_heuristic(data, rt_s, e_s, S, pi_s; rp_common...)
    end
    RP = actual_profit_of(rp_res)

    # ----- WS : wait-and-see (each scenario solved with a free first stage) ---
    # Solve each scenario alone (S = {sc}); its first stage adapts perfectly to
    # that scenario. Probability-weight the per-scenario actual profits.
    WS = 0.0
    wsum = sum(pi_s[sc] for sc in S)
    for sc in S
        Random.seed!(BASE_SEED + sc)
        ws_res = quiet() do
            stochastic_heuristic(data, rt_s, e_s, [sc], Dict(sc => 1.0); common...)
        end
        WS += (pi_s[sc]/wsum) * actual_profit_of(ws_res)
    end

    # ----- EEV : expected result of the expected-value solution --------------
    # 1) Solve the mean (EV) scenario to get the EV first-stage decision.
    #    Pool B is suppressed here (demand_rho = 0) so the EV commitment is made
    #    on the pre-booked Pool A demand under average weather only — Pool B is
    #    opportunistic on-demand traffic a mean-weather planner would not commit
    #    a fleet around. Pool B still appears in the per-scenario evaluation.
    rt_m, e_m, S_m, pi_m = build_mean_scenario(rt_s, e_s, S, pi_s)
    ev_common = merge(common, (demand_rho = 0.0,))
    Random.seed!(BASE_SEED + 999)
    ev_res = quiet() do
        stochastic_heuristic(data, rt_m, e_m, S_m, pi_m; ev_common...)
    end
    # 2) Fix that first-stage decision and evaluate it across the REAL scenarios.
    ev_fsd = FirstStageDecision(ev_res.accepted_passengers,
                                ev_res.active_evtols,
                                ev_res.tentative_sol,
                                ev_res.first_stage_cost)
    # Evaluate across real scenarios: probability-weighted true-fare profit, with
    # opening cost charged once.
    exp2nd = 0.0
    for sc in S
        sc_data = rp_res.scenario_results[sc].sc_data          # real scenario data
        rt_mat  = let Vmax = maximum(sc_data.V)
            m = zeros(Int, Vmax, Vmax)
            for (ij,v) in sc_data.rt; m[ij[1],ij[2]] = Int(round(v)); end
            m
        end
        r = quiet() do
            solve_second_stage(ev_fsd, sc_data, rt_mat;
                               maxTurnaround=MAXTURN, MaxTime_2nd=p.MaxTime_2nd_final,
                               top_c=TOP_C, price_boost=PRICE_BOOST,
                               hard_penalty=HARD_PENALTY)
        end
        rev = 0.0
        for ass in r.assignments
            a = ass.group; i = sc_data.op[a]; j = sc_data.dp[a]
            rev += sc_data.q[a]*(sc_data.fd[(i,j)]*(1-sc_data.so[a]) +
                                 sc_data.fs[(i,j)]*sc_data.so[a])
        end
        opcost = 0.0
        for plane in r.sol.planes
            for k in 1:plane.flightLegs
                opcost += sc_data.c[(Int(plane.route[k]), Int(plane.route[k+1]))]
            end
        end
        exp2nd += (pi_s[sc]/wsum) * (rev - opcost)
    end
    EEV = exp2nd - ev_res.first_stage_cost

    EVPI = WS - RP        # value of perfect information   (>= 0 expected)
    VSS  = RP - EEV       # value of the stochastic soln   (>= 0 expected)

    row = (
        season=SEASON_NAME, scenario_id=inst.id,
        n_VP=inst.vp, n_eVTOL=inst.evtol, n_passengers=inst.pax,
        WS=WS, RP=RP, EEV=EEV, EVPI=EVPI, VSS=VSS,
        ordering_ok = (WS + 1e-6 >= RP) && (RP + 1e-6 >= EEV),
    )
    push!(rows, row)
    @printf("  WS=%+.2f  RP=%+.2f  EEV=%+.2f  |  EVPI=%.2f  VSS=%.2f  %s\n",
            WS, RP, EEV, EVPI, VSS,
            row.ordering_ok ? "" : "  <-- ordering violated, check!")
    # When resuming (START_FROM > 1) append to the existing CSV so the already-
    # computed earlier instances are preserved; the header is written only if the
    # file does not yet exist. On a fresh run (START_FROM == 1) overwrite as usual.
    if START_FROM > 1 && isfile(out_path)
        CSV.write(out_path, DataFrame([row]); append = true)
    elseif START_FROM > 1
        CSV.write(out_path, DataFrame([row]))   # file gone — start it with header
    else
        CSV.write(out_path, DataFrame(rows))    # fresh run: overwrite
    end
end

println("\n" * "="^72)
println("Done. $(length(rows)) rows written to:\n  $out_path")
println("="^72)
