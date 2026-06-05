###############################################################################
# compute_vss_evpi.jl
#
# Computes the three standard stochastic-programming measures using the existing
# heuristic machinery:
#
#   RP   (Recourse Problem)            = E[profit] of the stochastic solution
#                                        (the here-and-now decision that hedges
#                                         across all scenarios). This is what
#                                         stochastic_heuristic already returns.
#
#   EEV  (Expected result of the       = take the first-stage decision that is
#         Expected Value solution)       optimal for the MEAN scenario, fix it,
#                                        and evaluate it across all scenarios.
#
#   WS   (Wait-and-See)                = for each scenario, solve it as if its
#                                        weather were known in advance; take the
#                                        probability-weighted average.
#
#   VSS  = RP  - EEV     (value of using stochastic programming)
#   EVPI = WS  - RP      (value of a perfect forecast)
#
# Required ordering (sanity check):  EEV <= RP <= WS
#
# This script assumes the heuristic source files and data are already loaded,
# i.e. run it AFTER MainCall_stochastic.jl has defined:
#   data, rt_s, e_s, S, pi_s, scenario_cache, rt_mats,
#   maxTurnaround, MaxTime_1st, MaxTime_2nd_search, MaxTime_2nd_final,
#   top_c, price_boost, hard_penalty
#
# If those are not in scope, set them here or include this from MainCall.
###############################################################################

using Printf

"""
    mean_scenario_index(S, pi_s)

Return the single scenario whose conditions are closest to the probability-
weighted "average" weather. We approximate the expected-value problem by the
most-probable scenario; if you have an explicit mean/nominal scenario, use that
instead. (A cleaner EEV uses a constructed mean-rt/mean-e scenario — see
make_mean_scenario below for that option.)
"""
function most_probable_scenario(S, pi_s)
    best_sc = first(S); best_p = -Inf
    for sc in S
        if pi_s[sc] > best_p
            best_p = pi_s[sc]; best_sc = sc
        end
    end
    return best_sc
end

"""
    make_mean_scenario(data, rt_s, e_s, S, pi_s)

Build the "mean" scenario by taking probability-weighted averages of rt and e
over all scenarios, then reusing make_scenario_data with a synthetic scenario
index (0) so the resulting sc_data has exactly the same shape as the real
scenario cache entries. Returns (mean_sc_data, mean_rt_mat).
"""
function make_mean_scenario(data, rt_s, e_s, S, pi_s)
    rt_aug = Dict{Tuple{Int,Int,Int},Float64}()
    e_aug  = Dict{Tuple{Int,Int,Int},Float64}()
    # rt_s / e_s only contain off-diagonal (i != j) keys — make_scenario_data
    # handles the i == j diagonal itself. Mirror that here: only average the
    # off-diagonal pairs, otherwise rt_s[(sc,i,i)] is a missing key.
    for i in data.V, j in data.V
        i == j && continue
        rt_avg = 0.0; e_avg = 0.0
        for sc in S
            rt_avg += pi_s[sc] * rt_s[(sc,i,j)]
            e_avg  += pi_s[sc] * e_s[(sc,i,j)]
        end
        rt_aug[(0,i,j)] = rt_avg
        e_aug[(0,i,j)]  = e_avg
    end
    mean_sc_data = make_scenario_data(data, 0, rt_aug, e_aug)

    Vmax = maximum(mean_sc_data.V)
    mean_rt_mat = zeros(Int, Vmax, Vmax)
    for (ij, v) in mean_sc_data.rt
        mean_rt_mat[ij[1], ij[2]] = v
    end
    return mean_sc_data, mean_rt_mat
end

"""
    compute_ws(scenario_cache, rt_mats, S, pi_s; kwargs...)

Wait-and-See: for each scenario, build the BEST first-stage decision tailored to
that scenario alone (perfect foresight), evaluate its profit in that scenario,
and take the probability-weighted average.

We obtain the per-scenario best decision by running the deterministic solve on
that scenario's data, extracting its first-stage decision, then scoring it in
that scenario only.
"""
function compute_ws(data, scenario_cache, rt_mats, S, pi_s;
                    maxTurnaround, MaxTime_det, top_c, price_boost, hard_penalty)
    ws = 0.0
    per_scenario = Dict{Int,Float64}()
    for sc in S
        sc_data = scenario_cache[sc]
        sc_rt   = rt_mats[sc]

        # Solve this scenario as a deterministic problem (perfect foresight).
        _, tent_sol, _ = HeuristicSA(maxTurnaround, MaxTime_det, sc_data, sc_rt, top_c)

        # The wait-and-see value of a scenario is simply the best achievable
        # profit in that scenario. true_det_objective scores the assembled plan
        # with real fares, costs, opening cost and feasibility penalties.
        prof, _ = true_det_objective(tent_sol, sc_data, sc_rt)

        per_scenario[sc] = prof
        ws += pi_s[sc] * prof
    end
    return ws, per_scenario
end

"""
    compute_eev(data, mean_rt, mean_e, scenario_cache, rt_mats, S, pi_s; ...)

Expected result of the Expected-Value solution:
  1. Solve the deterministic problem on the MEAN scenario.
  2. Extract its first-stage decision (accepted passengers + active eVTOLs).
  3. FIX that decision and evaluate it across ALL scenarios (recourse allowed).
"""
function compute_eev(data, mean_rt_mat, mean_sc_data,
                     scenario_cache, rt_mats, S, pi_s;
                     maxTurnaround, MaxTime_det, MaxTime_2nd,
                     top_c, price_boost, hard_penalty)
    # 1. Deterministic solve on the mean scenario.
    _, mean_sol, _ = HeuristicSA(maxTurnaround, MaxTime_det, mean_sc_data, mean_rt_mat, top_c)

    # 2. Extract the first-stage decision from the mean solution.
    #    Pass rt_mats so the robustness filter is applied across all scenarios,
    #    exactly as the stochastic heuristic does.
    eev_fsd = extract_first_stage_decision(mean_sol, data, mean_rt_mat, rt_mats)

    # 3. Evaluate that fixed decision across all scenarios.
    eev, _ = expected_objective(
        eev_fsd, scenario_cache, rt_mats, S, pi_s;
        maxTurnaround = maxTurnaround,
        MaxTime_2nd   = MaxTime_2nd,
        top_c         = top_c,
        price_boost   = price_boost,
        hard_penalty  = hard_penalty,
        full_output   = false,
    )
    return eev, eev_fsd
end

"""
    report_vss_evpi(RP, EEV, WS)

Print the three measures and the derived VSS and EVPI with the ordering check.
"""
function report_vss_evpi(RP, EEV, WS)
    VSS  = RP - EEV
    EVPI = WS - RP
    println("\n" * "="^72)
    println("STOCHASTIC PROGRAMMING MEASURES")
    println("="^72)
    @printf("  EEV  (expected-value solution, eval'd stochastically) : %+12.2f\n", EEV)
    @printf("  RP   (recourse / stochastic solution)                 : %+12.2f\n", RP)
    @printf("  WS   (wait-and-see, perfect foresight)                : %+12.2f\n", WS)
    println("  " * "-"^68)
    @printf("  VSS  = RP - EEV  (value of stochastic solution)       : %+12.2f\n", VSS)
    @printf("  EVPI = WS - RP   (expected value of perfect info)     : %+12.2f\n", EVPI)
    println("  " * "-"^68)
    ok = (EEV <= RP + 1e-6) && (RP <= WS + 1e-6)
    if ok
        println("  ✓ Ordering EEV <= RP <= WS holds.")
    else
        println("  ⚠ Ordering EEV <= RP <= WS VIOLATED — check search budgets/noise.")
        println("    (Heuristic noise can cause small violations; large ones signal a bug.)")
    end
    println("="^72)
    return (EEV = EEV, RP = RP, WS = WS, VSS = VSS, EVPI = EVPI)
end

###############################################################################
# Driver
#
# Call this after MainCall_stochastic.jl has run the stochastic heuristic, so
# you can pass in RP (the E[profit] it reported). It computes EEV and WS using
# the same data/parameters and prints VSS and EVPI.
#
# Use a GENEROUS deterministic budget (MaxTime_det) for WS and EEV so the
# per-scenario optima are well approximated — otherwise heuristic noise can make
# WS < RP, violating the ordering. WS especially should be solved well, since it
# is meant to represent the best achievable per scenario.
###############################################################################

function run_vss_evpi(data, rt_s, e_s, S, pi_s, scenario_cache, rt_mats, RP;
                      maxTurnaround   = 100,
                      MaxTime_det     = Int32(60),   # per-scenario deterministic budget
                      MaxTime_2nd     = Int32(30),   # recourse budget when evaluating EEV
                      top_c           = 4,
                      price_boost     = 8.0,
                      hard_penalty    = 50_000.0)

    println("\n[VSS/EVPI] Building mean scenario …")
    mean_sc_data, mean_rt_mat = make_mean_scenario(data, rt_s, e_s, S, pi_s)

    println("[VSS/EVPI] Computing EEV (expected-value solution, evaluated stochastically) …")
    EEV, eev_fsd = compute_eev(
        data, mean_rt_mat, mean_sc_data, scenario_cache, rt_mats, S, pi_s;
        maxTurnaround = maxTurnaround, MaxTime_det = MaxTime_det,
        MaxTime_2nd = MaxTime_2nd, top_c = top_c,
        price_boost = price_boost, hard_penalty = hard_penalty)
    println("           EEV first-stage accepted = $(sort(collect(eev_fsd.accepted_passengers)))")
    println("           EEV first-stage active   = $(sort(collect(eev_fsd.active_evtols)))")

    println("[VSS/EVPI] Computing WS (wait-and-see, perfect foresight per scenario) …")
    WS, ws_per = compute_ws(
        data, scenario_cache, rt_mats, S, pi_s;
        maxTurnaround = maxTurnaround, MaxTime_det = MaxTime_det,
        top_c = top_c, price_boost = price_boost, hard_penalty = hard_penalty)

    return report_vss_evpi(RP, EEV, WS)
end
