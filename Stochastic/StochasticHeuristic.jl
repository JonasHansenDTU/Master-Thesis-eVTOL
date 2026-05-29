###############################################################################
# StochasticHeuristic.jl
#
# Two-stage stochastic heuristic for eVTOL routing.
#
# ═══════════════════════════════════════════════════════════════════════════
# ARCHITECTURE
# ═══════════════════════════════════════════════════════════════════════════
#
# FIRST-STAGE DECISIONS (scenario-independent, made before weather is known)
#   1. accepted_passengers :: Set{Int}   — which passenger groups are committed
#   2. active_evtols       :: Set{Int}   — which eVTOLs are activated
#
# These are the ONLY binding first-stage outputs.
# A tentative plan is built to derive them, but is NOT binding.
#
# FIRST-STAGE COSTS (paid once, before any scenario is realized)
#   - Opening cost for each activated eVTOL
#   (Revenue is NOT counted first-stage — it is only earned when service
#    is actually delivered in the second stage)
#
# SECOND-STAGE DECISIONS (scenario-dependent, after weather is realized)
#   Given: accepted_passengers, active_evtols, scenario s
#   Decide: routes, timing, aircraft-to-passenger assignment
#   Hard constraint: ALL accepted passengers must be served
#   Optional: serve additional non-accepted passengers if capacity allows
#
# EXPECTED OBJECTIVE
#   E[profit] = - opening_cost * |active_evtols|
#               + Σ_s π_s * second_stage_profit(accepted, active_evtols, s)
#
#   where second_stage_profit = revenue_from_served - operating_costs
#                               - hard_penalty * unserved_accepted
#
# NON-ANTICIPATIVITY
#   accepted_passengers and active_evtols carry no scenario index.
#   All routing, timing, and assignment decisions are scenario-dependent.
#
# ═══════════════════════════════════════════════════════════════════════════


###############################################################################
# Helper: derive time_per_km from loaded data
###############################################################################

function derive_time_per_km(data)
    for (ij, d) in data.dist
        if d > 1.0 && haskey(data.rt, ij)
            return Float64(data.rt[ij]) / d
        end
    end
    error("Could not derive time_per_km: no non-zero-distance route pair found.")
end


###############################################################################
# First-stage decision object
#
# FirstStageDecision is the binding output of the first stage.
# It contains only the scenario-independent commitments.
# The tentative_plan is kept for reference and as a warm-start hint,
# but is never passed to the second stage as a constraint.
###############################################################################

struct FirstStageDecision
    accepted_passengers :: Set{Int}       # committed passenger groups
    active_evtols       :: Set{Int}       # activated eVTOL indices
    tentative_plan      :: allPlaneSolution  # reference only — not binding
    first_stage_cost    :: Float64        # opening costs already paid
end

"""
    extract_first_stage_decision(sol, data, rt_mat) -> FirstStageDecision

Given a tentative allPlaneSolution, extract:
  - which passengers are served (→ accepted)
  - which eVTOLs fly at least one leg (→ active)
  - the total opening cost already committed

The tentative plan itself is stored for reference but is NOT binding
in the second stage.
"""
function extract_first_stage_decision(sol::allPlaneSolution, data,
                                       rt_mat::Matrix{Int},
                                       rt_mats_all::Dict{Int,Matrix{Int}} = Dict{Int,Matrix{Int}}())

    assignments, _ = assign_passengersV2(sol, data, rt_mat)
    accepted_raw = Set{Int}(ass.group for ass in assignments)

    # Robustness filter: only accept a passenger if they can be served
    # feasibly across ALL scenarios.
    # Check 1: dt[a] + worst-case rt(op→dp) must fit before ET
    # Check 2: the eVTOL must also be able to return to an end vertiport
    #          after the trip: dt[a] + rt(op→dp) + rt(dp→end) <= ET
    # Check 3: departure time must be reachable from a base vertiport
    # Check 4: departure window must close before ET with margin for te
    #          (prevents very late passengers where dt[a] + w > ET)
    if !isempty(rt_mats_all)
        accepted = Set{Int}()
        for a in accepted_raw
            i = data.op[a]; j = data.dp[a]

            # Check 4 first (cheap, scenario-independent)
            # Passenger must depart within operating day including wait window
            if data.dt[a] > data.ET - Int(round(data.te))
                continue
            end

            robust = true
            for sc in keys(rt_mats_all)
                rt_mat_sc = rt_mats_all[sc]
                trip_rt   = rt_mat_sc[i, j]

                # Check 1: trip fits before ET
                if data.dt[a] + trip_rt > data.ET
                    robust = false; break
                end

                # Check 2: after serving the passenger (arriving at dp[a]=j),
                # at least one eVTOL must be able to reach one of its allowed
                # ending vertiports (end_vp[n]) before ET.
                can_return = any(
                    data.dt[a] + trip_rt +
                    minimum(rt_mat_sc[j, ev] for ev in data.end_vp[n];
                            init = typemax(Int)) <= data.ET
                    for n in data.N
                )
                if !can_return
                    robust = false; break
                end

                # Check 3: some eVTOL can reach op[a] before dt[a]
                reachable = any(rt_mat_sc[data.bv[n], i] <= data.dt[a]
                                for n in data.N)
                if !reachable
                    robust = false; break
                end
            end

            robust && push!(accepted, a)
        end
    else
        accepted = accepted_raw
    end

    active = Set{Int}(
        n for (n, plane) in enumerate(sol.planes)
        if plane.flightLegs > 0
    )

    first_stage_cost = length(active) * data.opening_cost

    return FirstStageDecision(accepted, active, sol, first_stage_cost)
end


###############################################################################
# Scenario-specific data object
###############################################################################

function make_scenario_data(data, sc::Int, rt_s, e_s)
    sc_rt   = Dict{Tuple{Int,Int}, Int}()
    sc_e    = Dict{Tuple{Int,Int}, Float64}()
    sc_dist = Dict{Tuple{Int,Int}, Float64}()

    for i in data.V, j in data.V
        if i == j
            sc_rt[(i,i)]   = 0
            sc_e[(i,i)]    = 0.0
            sc_dist[(i,i)] = 0.0
        else
            sc_rt[(i,j)]   = Int(ceil(rt_s[(sc,i,j)]))
            sc_e[(i,j)]    = e_s[(sc,i,j)]
            sc_dist[(i,j)] = e_s[(sc,i,j)] / data.battery_per_km
        end
    end

    return (
        infra          = data.infra,
        pax            = data.pax,
        plane          = data.plane,
        V              = data.V,
        A              = data.A,
        N              = data.N,
        M              = data.M,
        M_no0          = data.M_no0,
        M_mid          = data.M_mid,
        M_no_last      = data.M_no_last,
        T              = data.T,
        T_no0          = data.T_no0,
        bv             = data.bv,
        lat            = data.lat,
        lon            = data.lon,
        dist           = sc_dist,
        fd             = data.fd,
        fs             = data.fs,
        c              = data.c,
        e              = sc_e,
        rt             = sc_rt,
        end_vp         = data.end_vp,
        op             = data.op,
        dp             = data.dp,
        dt             = data.dt,
        q              = data.q,
        so             = data.so,
        p              = data.p,
        d              = data.d,
        cap_v          = data.cap_v,
        cap_u          = data.cap_u,
        opening_cost   = data.opening_cost,
        bmax           = data.bmax,
        bmid           = data.bmid,
        b_penalty      = data.b_penalty,
        bmin           = data.bmin,
        ec             = data.ec,
        te             = data.te,
        w              = data.w,
        ET             = data.ET,
        M1             = data.M1,
        M2a            = data.M2a,
        M2b            = data.M2b,
        M2c            = data.M2c,
        M3             = data.M3,
        battery_per_km = data.battery_per_km,
    )
end


###############################################################################
# Scenario data cache
###############################################################################

function build_scenario_cache(data, rt_s, e_s, S)
    cache   = Dict{Int, NamedTuple}()
    rt_mats = Dict{Int, Matrix{Int}}()
    for sc in S
        sc_data     = make_scenario_data(data, sc, rt_s, e_s)
        cache[sc]   = sc_data
        Vmax        = maximum(sc_data.V)
        mat         = zeros(Int, Vmax, Vmax)
        for (ij, v) in sc_data.rt
            mat[ij[1], ij[2]] = v
        end
        rt_mats[sc] = mat
    end
    return cache, rt_mats
end


###############################################################################
# Second-stage priority data
#
# Builds a modified scenario data object for the second-stage heuristic:
#
# 1. Only active eVTOLs are available (N restricted to active_evtols).
#    This enforces the first-stage fleet commitment: the second stage
#    cannot activate new aircraft.
#
# 2. fd/fs boosted for accepted passenger OD pairs so the construction
#    heuristic builds routes through them first.
#
# 3. p[a] set to hard_penalty for accepted passengers, so the fitness
#    function heavily penalises leaving them unserved.
#
# 4. opening_cost set to 0 in second-stage data because fleet activation
#    costs are already paid in the first stage — counting them again would
#    double-penalise.
###############################################################################

function make_second_stage_data(sc_data, fsd::FirstStageDecision;
                                  price_boost::Float64  = 50.0,
                                  hard_penalty::Float64 = 5_000.0)

    boosted_pairs = Set{Tuple{Int,Int}}(
        (sc_data.op[a], sc_data.dp[a]) for a in fsd.accepted_passengers
    )
    new_fd = copy(sc_data.fd)
    new_fs = copy(sc_data.fs)
    for (i, j) in boosted_pairs
        new_fd[(i,j)] = new_fd[(i,j)] * price_boost
        new_fs[(i,j)] = new_fs[(i,j)] * price_boost
    end

    new_p = copy(sc_data.p)
    for a in fsd.accepted_passengers
        new_p[a] = hard_penalty
    end

    return (
        infra          = sc_data.infra,
        pax            = sc_data.pax,
        plane          = sc_data.plane,
        V              = sc_data.V,
        A              = sc_data.A,
        N              = sc_data.N,        # ← all eVTOLs (activation cost paid in first stage)
        M              = sc_data.M,
        M_no0          = sc_data.M_no0,
        M_mid          = sc_data.M_mid,
        M_no_last      = sc_data.M_no_last,
        T              = sc_data.T,
        T_no0          = sc_data.T_no0,
        bv             = sc_data.bv,
        lat            = sc_data.lat,
        lon            = sc_data.lon,
        dist           = sc_data.dist,
        fd             = new_fd,
        fs             = new_fs,
        c              = sc_data.c,
        e              = sc_data.e,
        rt             = sc_data.rt,
        end_vp         = sc_data.end_vp,
        op             = sc_data.op,
        dp             = sc_data.dp,
        dt             = sc_data.dt,
        q              = sc_data.q,
        so             = sc_data.so,
        p              = new_p,
        d              = sc_data.d,
        cap_v          = sc_data.cap_v,
        cap_u          = sc_data.cap_u,
        opening_cost   = sc_data.opening_cost,  # keep real cost to discourage unnecessary activation
        bmax           = sc_data.bmax,
        bmid           = sc_data.bmid,
        b_penalty      = sc_data.b_penalty,
        bmin           = sc_data.bmin,
        ec             = sc_data.ec,
        te             = sc_data.te,
        w              = sc_data.w,
        ET             = sc_data.ET,
        M1             = sc_data.M1,
        M2a            = sc_data.M2a,
        M2b            = sc_data.M2b,
        M2c            = sc_data.M2c,
        M3             = sc_data.M3,
        battery_per_km = sc_data.battery_per_km,
    )
end


###############################################################################
# True second-stage profit
#
# Evaluated on ORIGINAL scenario data (not boosted) so numbers are realistic.
# Does NOT include opening costs (those are first-stage).
# Does NOT penalise non-accepted unserved passengers.
# DOES penalise unserved accepted passengers.
# DOES apply physical feasibility penalties (battery, timing, capacity).
###############################################################################

function second_stage_profit(sol::allPlaneSolution, sc_data,
                              sc_rt_mat::Matrix{Int},
                              accepted::Set{Int};
                              hard_penalty::Float64 = 5_000.0)

    assignments, _ = assign_passengersV2(sol, sc_data, sc_rt_mat)
    served = Set{Int}(ass.group for ass in assignments)

    value = 0.0

    # Revenue from served passengers (true prices)
    for ass in assignments
        a = ass.group
        i = sc_data.op[a]; j = sc_data.dp[a]
        value += sc_data.fd[(i,j)] * (1 - sc_data.so[a]) +
                 sc_data.fs[(i,j)] * sc_data.so[a]
    end

    # Operating cost for flown legs
    for plane in sol.planes
        for k in 1:plane.flightLegs
            value -= sc_data.c[(Int(plane.route[k]), Int(plane.route[k+1]))]
        end
    end

    # Physical feasibility penalties (battery, completion time, vertiport cap)
    P = FeasibilityCheck(
        Float32(sc_data.bmax), Float32(sc_data.bmid), Float32(sc_data.bmin),
        sc_data.dist, Float32(sc_data.ec), Float32(sc_data.battery_per_km),
        sol, sc_rt_mat,
        Int(round(sc_data.ET)), maximum(Int.(sc_data.T)),
        maximum(sc_data.V), sc_data.cap_v, sc_data.b_penalty
    )
    value -= sum(P)

    # Penalty only for unserved ACCEPTED passengers
    for a in accepted
        if !(a in served)
            value -= hard_penalty
        end
    end

    accepted_served = intersect(accepted, served)
    n_unserved      = length(accepted) - length(accepted_served)

    return value, assignments, accepted_served, n_unserved
end


###############################################################################
# Repair step: force-serve unserved accepted passengers
#
# After HeuristicSA, some accepted passengers may still be unserved because
# the heuristic did not find a route through their OD pair in time.
# This function scans every plane's route looking for legs that match
# the passenger's OD pair and adjusts turnaround times to fit their
# time window.  If a matching leg is found, the turnaround is adjusted
# so the passenger can be picked up.  If no matching leg exists at all,
# the passenger remains unserved — this is genuine infeasibility under
# this scenario and cannot be repaired without changing the route topology.
###############################################################################

function repair_unserved!(sol::allPlaneSolution, sc_data,
                           sc_rt_mat::Matrix{Int},
                           accepted::Set{Int})

    assignments, _ = assign_passengersV2(sol, sc_data, sc_rt_mat)
    served         = Set{Int}(ass.group for ass in assignments)
    still_unserved = Set{Int}()

    for a in accepted
        a in served && continue

        earliest = Int(round(sc_data.dt[a]))
        latest   = earliest + Int(round(sc_data.w))
        te       = Int(round(sc_data.te))
        fixed    = false

        for (pidx, plane) in enumerate(sol.planes)
            plane.flightLegs == 0 && continue

            cum_time = 0
            for k in 1:plane.flightLegs
                from = Int(plane.route[k])
                to   = Int(plane.route[k+1])

                if from == sc_data.op[a] && to == sc_data.dp[a]
                    # This leg has the right OD pair.
                    # Try to shift turnaround so departure falls in [earliest, latest].
                    min_ta  = k == 1 ? 0 : te
                    need_ta = earliest - cum_time
                    new_ta  = clamp(need_ta, min_ta, latest - cum_time)
                    new_dep = cum_time + new_ta

                    if earliest <= new_dep <= latest
                        plane.turnaroundTime[k] = Int32(new_ta)
                        fixed = true
                        break
                    end
                end

                cum_time += Int(plane.turnaroundTime[k]) +
                             sc_rt_mat[from, to]
            end
            fixed && break
        end

        fixed || push!(still_unserved, a)
    end

    return still_unserved
end


###############################################################################
# Second-stage solver
#
# Given a FirstStageDecision and a realized scenario:
# 1. Build second-stage data (restricted fleet, boosted prices, hard penalties)
# 2. Run HeuristicSA to find best operational plan
# 3. Run repair step to serve any missed accepted passengers
# 4. Evaluate true profit on original scenario data
###############################################################################

function solve_second_stage(fsd::FirstStageDecision,
                             sc_data, sc_rt_mat::Matrix{Int};
                             maxTurnaround::Int    = 100,
                             MaxTime_2nd::Int32    = Int32(10),
                             top_c::Int            = 4,
                             price_boost::Float64  = 50.0,
                             hard_penalty::Float64 = 5_000.0)

    ss_data = make_second_stage_data(sc_data, fsd;
                                      price_boost  = price_boost,
                                      hard_penalty = hard_penalty)

    # Run full heuristic — routes, times, aircraft assignment all free
    # (within the activated fleet)
    _, sc_sol, _ = HeuristicSA(maxTurnaround, MaxTime_2nd, ss_data,
                                sc_rt_mat, top_c)

    # Repair step: try to force-serve any accepted passengers still missed
    # by adjusting turnaround times on existing legs with the right OD pair.
    repair_unserved!(sc_sol, sc_data, sc_rt_mat, fsd.accepted_passengers)

    # Evaluate true profit on original scenario data
    profit, assignments, accepted_served, n_unserved = second_stage_profit(
        sc_sol, sc_data, sc_rt_mat, fsd.accepted_passengers;
        hard_penalty = hard_penalty
    )

    return (
        sol             = sc_sol,
        profit          = profit,
        assignments     = assignments,
        accepted_served = accepted_served,
        n_unserved      = n_unserved,
        sc_data         = sc_data,
    )
end


###############################################################################
# Expected objective
#
# E[profit] = - first_stage_cost
#             + Σ_s π_s * second_stage_profit(fsd, scenario s)
#
# The first-stage cost (opening costs) is subtracted once.
# Second-stage profits are averaged over scenarios.
###############################################################################

function expected_objective(fsd::FirstStageDecision,
                             scenario_cache, rt_mats, S, pi_s;
                             maxTurnaround::Int    = 100,
                             MaxTime_2nd::Int32    = Int32(5),
                             top_c::Int            = 4,
                             price_boost::Float64  = 50.0,
                             hard_penalty::Float64 = 5_000.0,
                             full_output::Bool     = false)

    exp_second_stage = 0.0
    scenario_results = full_output ? Dict{Int, NamedTuple}() : nothing

    for sc in S
        result = solve_second_stage(
            fsd, scenario_cache[sc], rt_mats[sc];
            maxTurnaround = maxTurnaround,
            MaxTime_2nd   = MaxTime_2nd,
            top_c         = top_c,
            price_boost   = price_boost,
            hard_penalty  = hard_penalty,
        )
        exp_second_stage += pi_s[sc] * result.profit
        if full_output
            scenario_results[sc] = result
        end
    end

    # First-stage cost is paid once (not per scenario)
    total = -fsd.first_stage_cost + exp_second_stage

    return total, scenario_results
end


###############################################################################
# True deterministic objective (for first-stage reporting)
###############################################################################

function true_det_objective(sol::allPlaneSolution, data, rt_mat::Matrix{Int})
    assignments, _ = assign_passengersV2(sol, data, rt_mat)
    value = 0.0

    for ass in assignments
        a = ass.group
        i = data.op[a]; j = data.dp[a]
        value += data.fd[(i,j)] * (1 - data.so[a]) + data.fs[(i,j)] * data.so[a]
    end
    for plane in sol.planes
        for k in 1:plane.flightLegs
            value -= data.c[(Int(plane.route[k]), Int(plane.route[k+1]))]
        end
        if plane.flightLegs > 0
            value -= data.opening_cost
        end
    end
    P = FeasibilityCheck(
        Float32(data.bmax), Float32(data.bmid), Float32(data.bmin),
        data.dist, Float32(data.ec), Float32(data.battery_per_km),
        sol, rt_mat,
        Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V),
        data.cap_v, data.b_penalty
    )
    value -= sum(P)

    return value, assignments
end


###############################################################################
# Main stochastic heuristic
#
# SEARCH LOOP
#   For each restart:
#     1. Run HeuristicSA on deterministic data → tentative plan
#     2. Extract FirstStageDecision (accepted passengers + active eVTOLs)
#     3. Evaluate E[profit] = -opening_costs + E[second_stage_profit]
#     4. Keep the FirstStageDecision with highest E[profit]
#
#   Final pass: re-solve all scenarios with longer budget for clean reporting.
###############################################################################

function stochastic_heuristic(data, rt_s, e_s, S, pi_s;
                               maxTurnaround::Int        = 100,
                               MaxTime_1st::Int32        = Int32(30),
                               MaxTime_2nd_search::Int32 = Int32(5),
                               MaxTime_2nd_final::Int32  = Int32(15),
                               top_c::Int                = 4,
                               price_boost::Float64      = 50.0,
                               hard_penalty::Float64     = 5_000.0,
                               n_restarts::Int           = 3)

    println("\n" * "="^72)
    println("TWO-STAGE STOCHASTIC HEURISTIC")
    println("="^72)
    println("  Scenarios              : $(length(S))")
    println("  First-stage budget     : $(MaxTime_1st)s × $n_restarts restarts")
    println("  Second-stage budget    : $(MaxTime_2nd_search)s (search), " *
            "$(MaxTime_2nd_final)s (final)")
    println("  Price boost            : $(price_boost)×")
    println("  Hard penalty/unserved  : $(hard_penalty)")
    println()
    println("  First-stage decisions  : accepted passengers + active eVTOLs")
    println("  First-stage cost       : opening_cost × |active_evtols| (paid once)")
    println("  Second-stage decisions : routes, timing, assignment (per scenario)")
    println("="^72)

    # Deterministic rt matrix for first-stage heuristic
    Vmax   = maximum(data.V)
    rt_det = zeros(Int, Vmax, Vmax)
    for i in data.V, j in data.V
        rt_det[i,j] = data.rt[(i,j)]
    end

    println("\n[Setup] Building scenario cache …")
    scenario_cache, rt_mats = build_scenario_cache(data, rt_s, e_s, S)
    println("[Setup] Done.")

    best_exp_obj = -Inf
    best_fsd     = nothing

    println("\n[First stage] Searching for best accepted-passenger + fleet commitment …")

    for restart in 1:n_restarts
        println("\n  ── Restart $restart / $n_restarts ──")

        # Build tentative plan under deterministic conditions
        _, tent_sol, _ = HeuristicSA(maxTurnaround, MaxTime_1st,
                                       data, rt_det, top_c)

        # Extract first-stage decision with robustness filter
        # (passes rt_mats so fragile passengers are excluded)
        fsd = extract_first_stage_decision(tent_sol, data, rt_det, rt_mats)

        tent_obj, _ = true_det_objective(tent_sol, data, rt_det)

        println("  Tentative det. obj     = $(round(tent_obj, digits=2))")
        println("  Accepted passengers    = $(sort(collect(fsd.accepted_passengers)))")
        println("  Active eVTOLs          = $(sort(collect(fsd.active_evtols)))")
        println("  First-stage cost       = $(round(fsd.first_stage_cost, digits=2))")

        # Evaluate E[profit] for this first-stage decision
        exp_obj, _ = expected_objective(
            fsd, scenario_cache, rt_mats, S, pi_s;
            maxTurnaround = maxTurnaround,
            MaxTime_2nd   = MaxTime_2nd_search,
            top_c         = top_c,
            price_boost   = price_boost,
            hard_penalty  = hard_penalty,
            full_output   = false,
        )

        println("  E[profit]              = $(round(exp_obj, digits=2))")

        if exp_obj > best_exp_obj
            best_exp_obj = exp_obj
            best_fsd     = fsd
            println("  *** New best E[profit] = $(round(best_exp_obj, digits=2))")
        end
    end

    println("\n[First stage] Best decision:")
    println("  Accepted : $(sort(collect(best_fsd.accepted_passengers)))")
    println("  Active   : $(sort(collect(best_fsd.active_evtols)))")
    println("  E[profit]: $(round(best_exp_obj, digits=2))")

    # Final second-stage pass with longer budget
    println("\n[Final] Re-solving all scenarios with $(MaxTime_2nd_final)s budget …")
    final_exp_obj, scenario_results = expected_objective(
        best_fsd, scenario_cache, rt_mats, S, pi_s;
        maxTurnaround = maxTurnaround,
        MaxTime_2nd   = MaxTime_2nd_final,
        top_c         = top_c,
        price_boost   = price_boost,
        hard_penalty  = hard_penalty,
        full_output   = true,
    )
    println("[Final] E[profit] = $(round(final_exp_obj, digits=2))")

    # Summary table
    println("\n" * "="^72)
    println("STOCHASTIC RESULTS SUMMARY")
    println("="^72)
    println("  Accepted passengers : $(sort(collect(best_fsd.accepted_passengers)))")
    println("  Active eVTOLs       : $(sort(collect(best_fsd.active_evtols)))")
    println()
    println("  Objective structure:")
    println(@sprintf("    First-stage cost (opening, paid once) : %+10.2f",
            -best_fsd.first_stage_cost))
    exp_ss = final_exp_obj + best_fsd.first_stage_cost
    println(@sprintf("    E[second-stage profit]                : %+10.2f", exp_ss))
    println(@sprintf("    E[total profit]                       : %+10.2f", final_exp_obj))
    println()
    println(@sprintf("  %-4s  %-18s  %8s  %12s  %10s",
            "sc", "label", "prob", "2nd_profit", "unserved"))
    println("  " * "-"^58)
    for sc in sort(collect(S))
        scen = SCENARIOS[sc]
        r    = scenario_results[sc]
        flag = r.n_unserved > 0 ? " !" : ""
        println(@sprintf("  %-4d  %-18s  %8.4f  %12.2f  %10d%s",
                sc, scen.label, pi_s[sc], r.profit, r.n_unserved, flag))
    end
    total_unserved = sum(scenario_results[sc].n_unserved for sc in S)
    println()
    if total_unserved == 0
        println("  ✓ All accepted passengers served in all scenarios.")
    else
        println("  WARNING: $total_unserved unserved accepted passengers across scenarios.")
        println("  Consider increasing price_boost or MaxTime_2nd_final.")
    end
    println()

    return (
        expected_obj        = final_exp_obj,
        accepted_passengers = best_fsd.accepted_passengers,
        active_evtols       = best_fsd.active_evtols,
        first_stage_cost    = best_fsd.first_stage_cost,
        tentative_sol       = best_fsd.tentative_plan,
        scenario_results    = scenario_results,
    )
end


###############################################################################
# Print results — compact version
###############################################################################

function print_stochastic_result(result, data)
    println("\n" * "="^72)
    println("FIRST-STAGE COMMITMENT")
    println("="^72)
    println("  Active eVTOLs: $(sort(collect(result.active_evtols)))")
    println(@sprintf("  Fleet activation cost: %.2f", result.first_stage_cost))
    println()
    println("  Accepted passengers:")
    println(@sprintf("  %-6s  %-8s  %-4s  %-8s  %-5s",
            "Group", "Route", "q", "dt", "so"))
    println("  " * "-"^38)
    for a in sort(collect(result.accepted_passengers))
        println(@sprintf("  %-6d  %d → %-4d  %-4d  %-8.0f  %-5d",
                a, data.op[a], data.dp[a], data.q[a], data.dt[a], data.so[a]))
    end

    println("\n" * "="^72)
    println("TENTATIVE FIRST-STAGE PLAN (reference only — not binding)")
    println("="^72)
    print_chromosome_table(result.tentative_sol)

    println("\n" * "="^72)
    println("SECOND-STAGE SUMMARY (per scenario)")
    println("="^72)
    println(@sprintf("  %-4s  %-18s  %8s  %12s  %10s  %s",
            "sc", "label", "prob", "2nd_profit", "unserved", "unserved groups"))
    println("  " * "-"^80)
    for sc in sort(collect(keys(result.scenario_results)))
        scen = SCENARIOS[sc]
        r    = result.scenario_results[sc]
        flag = r.n_unserved > 0 ? " !" : ""
        unserved_list = sort([a for a in result.accepted_passengers
                              if !(a in r.accepted_served)])
        unserved_str  = isempty(unserved_list) ? "—" : join(unserved_list, ", ")
        println(@sprintf("  %-4d  %-18s  %8.4f  %12.2f  %10d%s  [%s]",
                sc, scen.label, 1.0/length(keys(result.scenario_results)),
                r.profit, r.n_unserved, flag, unserved_str))
    end

    total_unserved = sum(result.scenario_results[sc].n_unserved
                         for sc in keys(result.scenario_results))
    println()
    if total_unserved == 0
        println("  ✓ All accepted passengers served in all scenarios.")
    else
        println("  WARNING: $total_unserved unserved accepted passengers across scenarios.")
    end

    println("\n" * "="^72)
    println("FINAL SUMMARY")
    println("="^72)
    @printf("  First-stage cost (opening)     = %+.2f\n", -result.first_stage_cost)
    exp_ss = result.expected_obj + result.first_stage_cost
    @printf("  E[second-stage profit]         = %+.2f\n", exp_ss)
    @printf("  E[total profit]  (E[profit])   = %+.2f\n", result.expected_obj)
    println("="^72)
end
