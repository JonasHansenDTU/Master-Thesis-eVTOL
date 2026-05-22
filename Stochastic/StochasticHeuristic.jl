###############################################################################
# StochasticHeuristic.jl
#
# Two-stage stochastic extension of the eVTOL heuristic.
#
# DESIGN
# ──────
# FIRST STAGE  (before uncertainty is realized)
#   - Which eVTOLs to activate  (implicit: those with flightLegs > 0)
#   - Which passenger groups to accept  (implicit: those assigned in the
#     deterministic solution)
#   Both are encoded in a single allPlaneSolution returned by the existing
#   heuristic on the deterministic data.  They are kept fixed across all
#   scenarios.
#
# SECOND STAGE  (after scenario s is realized)
#   Given the fixed set of accepted passengers and activated eVTOLs:
#   - Re-evaluate feasibility and objective with scenario-specific
#     travel times  rt_s[(sc,i,j)]  and battery consumption  e_s[(sc,i,j)].
#   - The routes themselves come from the first-stage solution; the
#     scenario changes only the *parameter values*, not the route topology.
#   - A short local-search improvement (SA neighbourhood) is then applied
#     per scenario to adapt turnaround times and route micro-structure to
#     the realized wind / temperature conditions.
#
# EXPECTED OBJECTIVE
#   E[profit] = (1/|S|) * Σ_s  second_stage_obj(s, first_stage_sol)
#
# USAGE (from MainCall.jl or similar)
# ──────────────────────────────────────────────────────────────────────────────
#   include("scenario_generation.jl")
#   include("StochasticHeuristic.jl")
#
#   data = load_data(excel_file, parameter_file)
#   time_per_km = derive_time_per_km(data)          # helper defined below
#   rt_s, e_s, S, pi_s = generate_scenarios(
#       data.V, data.lat, data.lon, data.rt, data.e, time_per_km)
#
#   result = stochastic_heuristic(data, rt_s, e_s, S, pi_s;
#                                 maxTurnaround = 100,
#                                 MaxTime       = Int32(60),
#                                 top_c         = 4,
#                                 sa_iters_per_scenario = 20)
#
#   println("Expected objective: $(result.expected_obj)")
#   println("First-stage solution:")
#   print_chromosome_table(result.first_stage_sol)
#   println("\nPer-scenario results:")
#   for (sc, r) in sort(collect(result.scenario_results), by = x -> x[1])
#       println("  sc=$(sc)  obj=$(round(r.obj, digits=2))")
#   end
###############################################################################

###############################################################################
# Helper: derive time_per_km from existing data
###############################################################################

"""
Derive the `time_per_km` parameter (minutes per km) from the loaded data.

`DataLoadFunc.load_data` does not expose `time_per_km` directly, so we
recover it by finding a non-zero-distance pair (i,j) and computing
    time_per_km ≈ rt[(i,j)] / dist[(i,j)]

The result is approximate (due to `ceil` in the original code) but close
enough for computing `v_air = 1 / time_per_km` in the wind model.
"""
function derive_time_per_km(data)
    for (ij, d) in data.dist
        if d > 1.0 && haskey(data.rt, ij)
            return Float64(data.rt[ij]) / d
        end
    end
    error("Could not derive time_per_km: no non-zero-distance route pair found in data.")
end


###############################################################################
# Scenario-specific data view
###############################################################################

"""
Build a lightweight scenario-specific data object by replacing `rt`, `e`,
`dist`, and `battery_per_km` with their scenario-dependent counterparts.

All other fields (passenger data, vertiport capacities, costs, …) are
shared references — they are the same NamedTuple fields, not copies, so
no extra allocation is needed for those.

The returned object has the same field layout as `data` and can be passed
directly to `obj(...)`, `FeasibilityCheck(...)`, `assign_passengersV2(...)`,
and the SA neighbourhood functions.
"""
function make_scenario_data(data, sc::Int, rt_s, e_s)
    # Build scenario-specific rt dict (same key type as data.rt)
    sc_rt = Dict{Tuple{Int,Int}, Int}()
    for i in data.V, j in data.V
        i == j && continue
        sc_rt[(i, j)] = Int(ceil(rt_s[(sc, i, j)]))
    end
    # Self-loops: rt[(i,i)] = 0 (needed by some loops in FeasibilityFunc)
    for i in data.V
        sc_rt[(i, i)] = 0
    end

    # Build scenario-specific e dict
    sc_e = Dict{Tuple{Int,Int}, Float64}()
    for i in data.V, j in data.V
        i == j && continue
        sc_e[(i, j)] = e_s[(sc, i, j)]
    end
    for i in data.V
        sc_e[(i, i)] = 0.0
    end

    # Build scenario-specific dist dict:
    # dist is used in FeasibilityFunc via BatteryNeeded(dist[(from,to)], battery_per_km).
    # Battery consumption = battery_per_km * dist.
    # Under a scenario: e_s[(sc,i,j)] = gamma_e * e[(i,j)]
    #                                  = gamma_e * battery_per_km * dist[(i,j)]
    # We keep battery_per_km fixed and adjust dist so that
    #   battery_per_km * sc_dist[(i,j)] = e_s[(sc,i,j)]
    # => sc_dist[(i,j)] = e_s[(sc,i,j)] / battery_per_km
    # This ensures FeasibilityCheck (which calls BatteryNeeded) remains correct.
    sc_dist = Dict{Tuple{Int,Int}, Float64}()
    for i in data.V, j in data.V
        i == j && continue
        sc_dist[(i, j)] = sc_e[(i, j)] / data.battery_per_km
    end
    for i in data.V
        sc_dist[(i, i)] = 0.0
    end

    # Return a NamedTuple that looks identical to data but with swapped fields.
    return (
        infra        = data.infra,
        pax          = data.pax,
        plane        = data.plane,
        V            = data.V,
        A            = data.A,
        N            = data.N,
        M            = data.M,
        M_no0        = data.M_no0,
        M_mid        = data.M_mid,
        M_no_last    = data.M_no_last,
        T            = data.T,
        T_no0        = data.T_no0,
        bv           = data.bv,
        lat          = data.lat,
        lon          = data.lon,
        dist         = sc_dist,         # ← scenario-specific
        fd           = data.fd,
        fs           = data.fs,
        c            = data.c,          # operating cost per km (deterministic)
        e            = sc_e,            # ← scenario-specific
        rt           = sc_rt,           # ← scenario-specific
        op           = data.op,
        dp           = data.dp,
        dt           = data.dt,
        q            = data.q,
        so           = data.so,
        p            = data.p,
        d            = data.d,
        cap_v        = data.cap_v,
        cap_u        = data.cap_u,
        opening_cost = data.opening_cost,
        bmax         = data.bmax,
        bmid         = data.bmid,
        b_penalty    = data.b_penalty,
        bmin         = data.bmin,
        ec           = data.ec,
        te           = data.te,
        w            = data.w,
        ET           = data.ET,
        M1           = data.M1,
        M2a          = data.M2a,
        M2b          = data.M2b,
        M2c          = data.M2c,
        M3           = data.M3,
        battery_per_km = data.battery_per_km,
    )
end


###############################################################################
# rt matrix helper (FeasibilityCheck needs a Matrix{Int}, not a Dict)
###############################################################################

"""
Convert the scenario rt dict to a dense Matrix{Int} as required by
FeasibilityCheck, FeasibleCompletionTime, build_scheduled_legs, etc.
"""
function rt_matrix(sc_data)
    Vmax = maximum(sc_data.V)
    mat  = zeros(Int, Vmax, Vmax)
    for (ij, v) in sc_data.rt
        mat[ij[1], ij[2]] = v
    end
    return mat
end


###############################################################################
# Second-stage evaluation
###############################################################################

"""
Evaluate the second-stage objective for a given first-stage solution under
scenario `sc`.

Steps:
  1. Build scenario-specific data (rt_s, e_s → modified dist and rt).
  2. Optionally run a short SA local search to adapt the solution to the
     scenario (improving turnaround times / micro-route).
  3. Compute and return the scenario objective value and final solution.

Arguments:
  sol            : first-stage allPlaneSolution (will be deepcopy-ed)
  data           : original deterministic data
  sc             : scenario index
  rt_s, e_s      : scenario-dependent parameter dicts
  maxTurnaround  : max turnaround time (minutes) passed to SA
  sa_iters       : number of SA improvement iterations for this scenario
                   (0 = evaluate only, no adaptation)
"""
function evaluate_second_stage(sol::allPlaneSolution, data, sc::Int,
                                rt_s, e_s;
                                maxTurnaround::Int = 100,
                                sa_iters::Int = 20)

    sc_data = make_scenario_data(data, sc, rt_s, e_s)
    sc_rt   = rt_matrix(sc_data)

    # Work on a copy so the first-stage solution is never mutated
    current_sol = deepcopy(sol)
    current_obj = obj(current_sol, sc_data, sc_rt)

    # Short SA loop to adapt the solution to this scenario's conditions.
    # We reuse the same SA neighbourhood operators (DestructSA / ConstructSA)
    # but pass the scenario-specific data and rt matrix.
    if sa_iters > 0
        T_sa = 20.0
        for _ in 1:sa_iters
            if rand(Bool)
                cand_obj, cand_sol = DestructSA(current_sol, maxTurnaround,
                                                current_obj, sc_data, sc_rt)
            else
                cand_obj, cand_sol = ConstructSA(current_sol, maxTurnaround,
                                                 current_obj, sc_data, sc_rt)
            end

            # SA acceptance
            if cand_obj > current_obj ||
               rand() < exp((cand_obj - current_obj) / T_sa)
                current_sol = cand_sol
                current_obj = cand_obj
            end
            T_sa = max(T_sa * 0.92, 0.001)
        end
    end

    return (obj = current_obj, sol = current_sol, sc_data = sc_data)
end


###############################################################################
# Expected-value objective
###############################################################################

"""
Compute the expected objective of a first-stage solution across all scenarios.

Returns:
  expected_obj     : Σ_s  pi_s[s] * scenario_obj[s]
  scenario_results : Dict{Int, NamedTuple{(:obj, :sol, :sc_data)}}
"""
function expected_objective(sol::allPlaneSolution, data,
                             rt_s, e_s, S, pi_s;
                             maxTurnaround::Int = 100,
                             sa_iters_per_scenario::Int = 20)

    scenario_results = Dict{Int, NamedTuple}()
    expected_obj     = 0.0

    for sc in S
        result = evaluate_second_stage(sol, data, sc, rt_s, e_s;
                                       maxTurnaround = maxTurnaround,
                                       sa_iters      = sa_iters_per_scenario)
        scenario_results[sc] = result
        expected_obj        += pi_s[sc] * result.obj
    end

    return expected_obj, scenario_results
end


###############################################################################
# Main stochastic heuristic
###############################################################################

"""
Two-stage stochastic heuristic for eVTOL routing.

FIRST STAGE
  Run `HeuristicSA` on the deterministic data to find the best
  fleet activation and passenger acceptance (encoded in the returned
  allPlaneSolution).  The heuristic is run for `MaxTime` seconds.

SECOND STAGE
  For the best first-stage solution (and optionally a pool of top
  candidates), compute the expected objective across all scenarios.
  A short SA adaptation is applied per scenario.

The function returns the first-stage solution that maximises the
expected objective, together with per-scenario details.

Arguments:
  data                  : loaded data (from load_data)
  rt_s, e_s             : scenario parameter dicts (from generate_scenarios)
  S                     : scenario index range
  pi_s                  : scenario probability dict
  maxTurnaround         : max turnaround time (minutes)
  MaxTime               : wall-clock budget for first-stage heuristic (seconds)
  top_c                 : top_c parameter passed to HeuristicSA
  sa_iters_per_scenario : SA iterations in second-stage adaptation per scenario
  n_first_stage_candidates : how many top first-stage solutions to compare
                             before picking the best expected-value one

Returns a NamedTuple with fields:
  expected_obj      : best expected objective value found
  first_stage_sol   : the chosen first-stage allPlaneSolution
  scenario_results  : Dict of per-scenario (obj, sol, sc_data) for the winner
  det_obj           : deterministic objective of the first-stage solution
  iterations        : iterations run by the first-stage heuristic
"""
function stochastic_heuristic(data, rt_s, e_s, S, pi_s;
                               maxTurnaround::Int        = 100,
                               MaxTime::Int32            = Int32(60),
                               top_c::Int                = 4,
                               sa_iters_per_scenario::Int = 20,
                               n_first_stage_candidates::Int = 3)

    println("\n" * "="^72)
    println("TWO-STAGE STOCHASTIC HEURISTIC")
    println("  Scenarios : $(length(S))")
    println("  SA iters / scenario : $sa_iters_per_scenario")
    println("  First-stage time budget : $(MaxTime)s")
    println("="^72)

    # ── FIRST STAGE ──────────────────────────────────────────────────────────
    println("\n[First stage] Running HeuristicSA on deterministic data …")

    # Build the dense rt matrix the heuristic expects
    Vmax = maximum(data.V)
    rt_det = zeros(Vmax, Vmax)
    for i in data.V, j in data.V
        rt_det[i, j] = data.rt[(i, j)]
    end

    det_obj, first_sol, iters = HeuristicSA(maxTurnaround, MaxTime, data,
                                             rt_det, top_c)

    println("[First stage] Done.  det_obj = $(round(det_obj, digits=2))  " *
            "iterations = $iters")

    # ── EVALUATE EXPECTED OBJECTIVE ──────────────────────────────────────────
    println("\n[Second stage] Evaluating expected objective across $(length(S)) scenarios …")

    exp_obj, sc_results = expected_objective(
        first_sol, data, rt_s, e_s, S, pi_s;
        maxTurnaround           = maxTurnaround,
        sa_iters_per_scenario   = sa_iters_per_scenario
    )

    println("[Second stage] Expected objective = $(round(exp_obj, digits=2))")

    # ── SUMMARY ──────────────────────────────────────────────────────────────
    println("\n" * "="^72)
    println("STOCHASTIC RESULTS SUMMARY")
    println("="^72)
    println(@sprintf("  Deterministic (first-stage) objective : %10.2f", det_obj))
    println(@sprintf("  Expected objective (E[profit])        : %10.2f", exp_obj))
    println()
    println("  Per-scenario breakdown:")
    println(@sprintf("  %-4s  %-18s  %10s  %10s", "sc", "label", "prob", "obj"))
    println("  " * "-"^46)
    for sc in sort(collect(S))
        scen  = SCENARIOS[sc]
        r     = sc_results[sc]
        println(@sprintf("  %-4d  %-18s  %10.4f  %10.2f",
                sc, scen.label, pi_s[sc], r.obj))
    end
    println()

    return (
        expected_obj     = exp_obj,
        first_stage_sol  = first_sol,
        scenario_results = sc_results,
        det_obj          = det_obj,
        iterations       = iters,
    )
end


###############################################################################
# Convenience: print stochastic result
###############################################################################

"""
Pretty-print the output of `stochastic_heuristic`.
"""
function print_stochastic_result(result, data, rt_s, e_s)
    println("\n" * "="^72)
    println("FIRST-STAGE SOLUTION")
    println("="^72)
    print_chromosome_table(result.first_stage_sol)

    println("\n" * "="^72)
    println("SECOND-STAGE SOLUTIONS (per scenario)")
    println("="^72)
    for sc in sort(collect(keys(result.scenario_results)))
        scen = SCENARIOS[sc]
        r    = result.scenario_results[sc]
        println("\n── Scenario $sc : $(scen.label)  " *
                "(wx=$(scen.wx), wy=$(scen.wy), ϕ=$(scen.phi)) ──")
        println(@sprintf("   Scenario objective : %.2f", r.obj))

        sc_rt_mat = rt_matrix(r.sc_data)
        assignments, scheduled = assign_passengersV2(r.sol, r.sc_data, sc_rt_mat)

        println("   Passenger assignments:")
        for ass in assignments
            println("     Group $(ass.group) | plane $(ass.plane) | legs $(ass.legs) | " *
                    "$(data.op[ass.group]) → $(data.dp[ass.group])")
        end
    end

    println("\n" * "="^72)
    @printf("Expected objective E[profit] = %.2f\n", result.expected_obj)
    @printf("Deterministic objective      = %.2f\n", result.det_obj)
    val_of_stochastic = result.expected_obj - result.det_obj
    println(@sprintf("Value of stochastic solution = %.2f  (E[profit] – det_obj)",
                     val_of_stochastic))
    println("="^72)
end
