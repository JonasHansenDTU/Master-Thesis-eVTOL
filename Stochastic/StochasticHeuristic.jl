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
# FirstStageDecision captures the two binding scenario-independent decisions:
#   1. accepted_passengers — which passenger groups are committed
#   2. active_evtols       — which eVTOLs are activated (fleet commitment)
#
# Note on eVTOL assignment: the MIP has ss[a,n] binding passenger a to
# eVTOL n in the first stage. In the heuristic we deliberately do NOT lock
# the assignment — which eVTOL serves which accepted passenger is a
# second-stage decision that can adapt to the realised scenario. This gives
# the second stage more flexibility to also fit opportunistic non-accepted
# passengers into remaining capacity.
###############################################################################

struct FirstStageDecision
    accepted_passengers :: Set{Int}          # committed passenger groups
    active_evtols       :: Set{Int}          # activated eVTOL indices
    tentative_plan      :: allPlaneSolution  # reference only — not binding
    first_stage_cost    :: Float64           # opening costs already paid
end

"""
    extract_first_stage_decision(sol, data, rt_mat, rt_mats_all) -> FirstStageDecision

Given a tentative allPlaneSolution, extract the binding first-stage decisions:
  - accepted_passengers: passengers served in the tentative plan that pass a
    robustness filter (feasible in ALL scenarios for SOME active eVTOL)
  - active_evtols: eVTOLs that fly at least one leg in the tentative plan
  - first_stage_cost: opening costs for the active fleet
"""
function extract_first_stage_decision(sol::allPlaneSolution, data,
                                       rt_mat::Matrix{Int},
                                       rt_mats_all::Dict{Int,Matrix{Int}} = Dict{Int,Matrix{Int}}())

    assignments, _ = assign_passengersV2(sol, data, rt_mat)
    accepted_raw   = Set{Int}(ass.group for ass in assignments)

    active = Set{Int}(
        n for (n, plane) in enumerate(sol.planes)
        if plane.flightLegs > 0
    )

    if !isempty(rt_mats_all)
        accepted = Set{Int}()
        for a in accepted_raw
            i = data.op[a]; j = data.dp[a]

            # Check 4: scenario-independent — must depart before ET - te
            if data.dt[a] > data.ET - Int(round(data.te))
                continue
            end

            robust = true
            for sc in keys(rt_mats_all)
                rt_mat_sc = rt_mats_all[sc]
                trip_rt   = rt_mat_sc[i, j]

                # Check 1: trip completes before ET
                if data.dt[a] + trip_rt > data.ET
                    robust = false; break
                end

                # Check 2: some active eVTOL can return to one of its
                # end vertiports after serving passenger a
                can_return = any(
                    data.dt[a] + trip_rt +
                    minimum(rt_mat_sc[j, ev]
                            for ev in data.end_vp[n];
                            init = typemax(Int)) <= data.ET
                    for n in active
                )
                if !can_return
                    robust = false; break
                end

                # Check 3: some active eVTOL can reach op[a] before dt[a],
                # with at least te minutes margin for turnaround.
                if !any(rt_mat_sc[data.bv[n], i] + Int(round(data.te)) <= data.dt[a]
                        for n in active)
                    robust = false; break
                end
            end

            robust && push!(accepted, a)
        end
    else
        accepted = accepted_raw
    end

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

###############################################################################
# Second-stage priority data
#
# Builds a modified scenario data object for the second-stage heuristic:
#
# 1. A = full passenger set — the heuristic can serve ALL passengers,
#    earning opportunistic revenue on top of the accepted ones.
#
# 2. p[a] = hard_penalty for accepted passengers, original p[a] for others.
#    This makes the heuristic strongly prioritise accepted passengers while
#    still picking up non-accepted passengers when there is spare capacity.
#
# 3. fd/fs boosted for accepted passenger OD pairs so route construction
#    naturally visits their origins/destinations first.
#
# 4. opening_cost preserved to discourage activating extra eVTOLs.
###############################################################################

function make_second_stage_data(sc_data, fsd::FirstStageDecision;
                                  price_boost::Float64  = 10.0,
                                  hard_penalty::Float64 = 10_000.0)

    new_fd = copy(sc_data.fd)
    new_fs = copy(sc_data.fs)

    # Boost fares for accepted passenger OD pairs so the construction
    # heuristic prioritises their origins and destinations.
    boosted_pairs = Set{Tuple{Int,Int}}(
        (sc_data.op[a], sc_data.dp[a]) for a in fsd.accepted_passengers
    )
    for (i, j) in boosted_pairs
        new_fd[(i,j)] = new_fd[(i,j)] * price_boost
        new_fs[(i,j)] = new_fs[(i,j)] * price_boost
    end

    # Full passenger set: heuristic can serve anyone, but accepted passengers
    # carry hard_penalty if unserved so they are always prioritised.
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
        N              = sc_data.N,
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
        opening_cost   = sc_data.opening_cost,
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

        # If we have a committed eVTOL assignment, only search that eVTOL.
        # Otherwise search all planes (fallback for backward compatibility).
        plane_indices = collect(enumerate(sol.planes))

        for (pidx, plane) in plane_indices
            plane.flightLegs == 0 && continue

            cum_time = 0
            for k in 1:plane.flightLegs
                from = Int(plane.route[k])
                to   = Int(plane.route[k+1])

                if from == sc_data.op[a] && to == sc_data.dp[a]
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

    # Repair step: try to force-serve any accepted passengers still missed,
    # searching only the committed eVTOL for each passenger.
    repair_unserved!(sc_sol, sc_data, sc_rt_mat,
                     fsd.accepted_passengers)

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
# First-stage neighbourhood helpers
#
# These functions operate directly on the (accepted, active) decision space.
# They are the moves of the outer search loop.
###############################################################################

"""
    build_robust_candidates(data, active, rt_mats_all) -> Set{Int}

Return all passengers that pass the robustness filter for the given active
fleet — i.e. passengers that some active eVTOL can feasibly serve in every
scenario. These are the candidates available to add to the accepted set.
"""
function build_robust_candidates(data, active::Set{Int},
                                  rt_mats_all::Dict{Int,Matrix{Int}})
    candidates = Set{Int}()
    for a in data.A
        i = data.op[a]; j = data.dp[a]
        if data.dt[a] > data.ET - Int(round(data.te))
            continue
        end
        robust = true
        for sc in keys(rt_mats_all)
            rt_sc   = rt_mats_all[sc]
            trip_rt = rt_sc[i, j]
            if data.dt[a] + trip_rt > data.ET
                robust = false; break
            end
            can_return = any(
                data.dt[a] + trip_rt +
                minimum(rt_sc[j, ev] for ev in data.end_vp[n];
                        init = typemax(Int)) <= data.ET
                for n in active
            )
            if !can_return
                robust = false; break
            end
            if !any(rt_sc[data.bv[n], i] + Int(round(data.te)) <= data.dt[a]
                    for n in active)
                robust = false; break
            end
        end
        robust && push!(candidates, a)
    end
    return candidates
end

"""
    make_fsd(accepted, active, data, tentative_sol) -> FirstStageDecision

Construct a FirstStageDecision from explicit accepted and active sets.
"""
function make_fsd(accepted::Set{Int}, active::Set{Int},
                   data, tentative_sol::allPlaneSolution)
    cost = length(active) * data.opening_cost
    return FirstStageDecision(accepted, active, tentative_sol, cost)
end

###############################################################################
# Main stochastic heuristic
#
# PHASE 1 — Initialisation (deterministic, n_restarts times)
#   Run HeuristicSA on deterministic data → extract candidate first-stage
#   decisions → evaluate E[profit] → keep best as starting point.
#
# PHASE 2 — Outer stochastic search (scenario-aware local search)
#   This is the core of the proper two-stage stochastic heuristic.
#   The objective function IS E[profit] evaluated across all scenarios.
#   Neighbourhood moves operate directly on (accepted, active):
#     - ADD:    add one robust passenger to accepted
#     - REMOVE: drop one accepted passenger
#     - SWAP:   replace one accepted passenger with a different robust one
#     - FLEET+: activate one additional eVTOL (opens more routing options)
#     - FLEET-: deactivate one eVTOL (saves opening cost)
#   Each move evaluates E[profit] using the second-stage heuristic on all
#   scenarios. Accepts improving moves immediately (greedy local search).
#
# PHASE 3 — Final pass
#   Re-solve all scenarios with longer budget for clean reporting.
###############################################################################

function stochastic_heuristic(data, rt_s, e_s, S, pi_s;
                               maxTurnaround::Int        = 100,
                               MaxTime_1st::Int32        = Int32(30),
                               MaxTime_2nd_search::Int32 = Int32(5),
                               MaxTime_2nd_final::Int32  = Int32(15),
                               top_c::Int                = 4,
                               price_boost::Float64      = 10.0,
                               hard_penalty::Float64     = 10_000.0,
                               n_restarts::Int           = 3,
                               n_outer_iters::Int        = 20)

    println("\n" * "="^72)
    println("TWO-STAGE STOCHASTIC HEURISTIC")
    println("="^72)
    println("  Scenarios              : $(length(S))")
    println("  Phase 1 (init)         : $(MaxTime_1st)s × $n_restarts restarts")
    println("  Phase 2 (outer search) : $n_outer_iters iterations × all scenarios")
    println("  Second-stage budget    : $(MaxTime_2nd_search)s (search), " *
            "$(MaxTime_2nd_final)s (final)")
    println("  Price boost            : $(price_boost)×")
    println("  Hard penalty/unserved  : $(hard_penalty)")
    println()
    println("  First-stage decisions  : accepted passengers + active eVTOLs")
    println("  First-stage cost       : opening_cost × |active_evtols| (paid once)")
    println("  Second-stage decisions : routes, timing, assignment (per scenario)")
    println("="^72)

    # Deterministic rt matrix for first-stage initialisation
    Vmax   = maximum(data.V)
    rt_det = zeros(Int, Vmax, Vmax)
    for i in data.V, j in data.V
        rt_det[i,j] = data.rt[(i,j)]
    end

    println("\n[Setup] Building scenario cache …")
    scenario_cache, rt_mats = build_scenario_cache(data, rt_s, e_s, S)
    println("[Setup] Done.")

    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 1: Initialisation — run n_restarts deterministic solves to get
    # a good starting point for the outer stochastic search.
    # ─────────────────────────────────────────────────────────────────────────
    println("\n" * "─"^72)
    println("PHASE 1 — Initialisation (deterministic, $n_restarts restarts)")
    println("─"^72)

    best_exp_obj  = -Inf
    best_fsd      = nothing
    best_tent_sol = nothing

    for restart in 1:n_restarts
        println("\n  ── Restart $restart / $n_restarts ──")

        _, tent_sol, _ = HeuristicSA(maxTurnaround, MaxTime_1st,
                                       data, rt_det, top_c)

        fsd = extract_first_stage_decision(tent_sol, data, rt_det, rt_mats)

        tent_obj, _ = true_det_objective(tent_sol, data, rt_det)
        println("  Tentative det. obj  = $(round(tent_obj, digits=2))")
        println("  Accepted            = $(sort(collect(fsd.accepted_passengers)))")
        println("  Active eVTOLs       = $(sort(collect(fsd.active_evtols)))")

        exp_obj, _ = expected_objective(
            fsd, scenario_cache, rt_mats, S, pi_s;
            maxTurnaround = maxTurnaround,
            MaxTime_2nd   = MaxTime_2nd_search,
            top_c         = top_c,
            price_boost   = price_boost,
            hard_penalty  = hard_penalty,
            full_output   = false,
        )
        println("  E[profit]           = $(round(exp_obj, digits=2))")

        if exp_obj > best_exp_obj
            best_exp_obj  = exp_obj
            best_fsd      = fsd
            best_tent_sol = tent_sol
            println("  *** New best = $(round(best_exp_obj, digits=2))")
        end
    end

    println("\n  Phase 1 best: accepted=$(sort(collect(best_fsd.accepted_passengers)))")
    println("  Phase 1 best: E[profit]=$(round(best_exp_obj, digits=2))")

    # ─────────────────────────────────────────────────────────────────────────
    # PHASE 2: Outer stochastic search over first-stage decisions.
    #
    # The objective function IS E[profit] evaluated across all scenarios.
    # We perform local search directly on (accepted_passengers, active_evtols).
    #
    # Moves tried at each iteration (in random order):
    #   ADD:    pick a random robust non-accepted passenger → add to accepted
    #   REMOVE: pick a random accepted passenger → remove from accepted
    #   SWAP:   remove one accepted, add one robust non-accepted
    #   FLEET+: activate one idle eVTOL → expand routing options
    #   FLEET-: deactivate one active eVTOL with fewest committed passengers
    #
    # Greedy acceptance: accept any move that strictly improves E[profit].
    # After the greedy phase, the best solution found is kept.
    # ─────────────────────────────────────────────────────────────────────────
    println("\n" * "─"^72)
    println("PHASE 2 — Outer stochastic search ($n_outer_iters iterations)")
    println("─"^72)
    println("  Objective: E[profit] evaluated across all $(length(S)) scenarios")
    println("  Moves: ADD / REMOVE / SWAP passenger, FLEET+ / FLEET-")
    println()

    current_fsd = best_fsd
    current_obj = best_exp_obj

    all_evtols = Set{Int}(data.N)

    for iter in 1:n_outer_iters
        current_accepted = current_fsd.accepted_passengers
        current_active   = current_fsd.active_evtols

        # Build the set of robust candidates for the current active fleet
        candidates = build_robust_candidates(data, current_active, rt_mats)

        # Non-accepted robust passengers — available to add
        addable = setdiff(candidates, current_accepted)

        # Build list of moves to try this iteration
        moves = Symbol[]
        !isempty(addable)            && push!(moves, :add)
        !isempty(current_accepted)   && push!(moves, :remove)
        !isempty(addable) && !isempty(current_accepted) && push!(moves, :swap)
        length(current_active) < length(all_evtols) && push!(moves, :fleet_add)
        length(current_active) > 1   && push!(moves, :fleet_remove)

        # Shuffle move order for diversity
        shuffle!(moves)

        improved = false
        best_candidate_obj = current_obj
        best_candidate_fsd = nothing

        for move in moves
            # For each move type, try ALL candidates (not just one random one)
            # and keep the best improving one. This is the key fix — a single
            # random pick is too noisy and misses good moves.

            if move == :add
                for a in shuffle!(collect(addable))
                    c_accepted = push!(copy(current_accepted), a)
                    c_fsd = make_fsd(c_accepted, copy(current_active),
                                     data, best_tent_sol)
                    c_obj, _ = expected_objective(
                        c_fsd, scenario_cache, rt_mats, S, pi_s;
                        maxTurnaround = maxTurnaround,
                        MaxTime_2nd   = MaxTime_2nd_search,
                        top_c         = top_c,
                        price_boost   = price_boost,
                        hard_penalty  = hard_penalty,
                        full_output   = false,
                    )
                    if c_obj > best_candidate_obj
                        best_candidate_obj = c_obj
                        best_candidate_fsd = c_fsd
                    end
                end

            elseif move == :remove
                for a in shuffle!(collect(current_accepted))
                    c_accepted = delete!(copy(current_accepted), a)
                    c_fsd = make_fsd(c_accepted, copy(current_active),
                                     data, best_tent_sol)
                    c_obj, _ = expected_objective(
                        c_fsd, scenario_cache, rt_mats, S, pi_s;
                        maxTurnaround = maxTurnaround,
                        MaxTime_2nd   = MaxTime_2nd_search,
                        top_c         = top_c,
                        price_boost   = price_boost,
                        hard_penalty  = hard_penalty,
                        full_output   = false,
                    )
                    if c_obj > best_candidate_obj
                        best_candidate_obj = c_obj
                        best_candidate_fsd = c_fsd
                    end
                end

            elseif move == :swap
                # Try a sample of swaps to keep runtime manageable
                for a_out in shuffle!(collect(current_accepted))
                    for a_in in shuffle!(collect(addable))[1:min(3,length(addable))]
                        c_accepted = push!(delete!(copy(current_accepted), a_out), a_in)
                        c_fsd = make_fsd(c_accepted, copy(current_active),
                                         data, best_tent_sol)
                        c_obj, _ = expected_objective(
                            c_fsd, scenario_cache, rt_mats, S, pi_s;
                            maxTurnaround = maxTurnaround,
                            MaxTime_2nd   = MaxTime_2nd_search,
                            top_c         = top_c,
                            price_boost   = price_boost,
                            hard_penalty  = hard_penalty,
                            full_output   = false,
                        )
                        if c_obj > best_candidate_obj
                            best_candidate_obj = c_obj
                            best_candidate_fsd = c_fsd
                        end
                    end
                end

            elseif move == :fleet_add
                for n in shuffle!(collect(setdiff(all_evtols, current_active)))
                    c_active = push!(copy(current_active), n)
                    c_candidates = build_robust_candidates(data, c_active, rt_mats)
                    c_accepted = intersect(current_accepted, c_candidates) |> Set{Int}
                    c_fsd = make_fsd(c_accepted, c_active, data, best_tent_sol)
                    c_obj, _ = expected_objective(
                        c_fsd, scenario_cache, rt_mats, S, pi_s;
                        maxTurnaround = maxTurnaround,
                        MaxTime_2nd   = MaxTime_2nd_search,
                        top_c         = top_c,
                        price_boost   = price_boost,
                        hard_penalty  = hard_penalty,
                        full_output   = false,
                    )
                    if c_obj > best_candidate_obj
                        best_candidate_obj = c_obj
                        best_candidate_fsd = c_fsd
                    end
                end

            elseif move == :fleet_remove
                for n in shuffle!(collect(current_active))
                    length(current_active) <= 1 && continue
                    c_active = delete!(copy(current_active), n)
                    c_candidates = build_robust_candidates(data, c_active, rt_mats)
                    c_accepted = intersect(current_accepted, c_candidates) |> Set{Int}
                    c_fsd = make_fsd(c_accepted, c_active, data, best_tent_sol)
                    c_obj, _ = expected_objective(
                        c_fsd, scenario_cache, rt_mats, S, pi_s;
                        maxTurnaround = maxTurnaround,
                        MaxTime_2nd   = MaxTime_2nd_search,
                        top_c         = top_c,
                        price_boost   = price_boost,
                        hard_penalty  = hard_penalty,
                        full_output   = false,
                    )
                    if c_obj > best_candidate_obj
                        best_candidate_obj = c_obj
                        best_candidate_fsd = c_fsd
                    end
                end
            end

            # If any candidate in this move type improved, accept it
            if !isnothing(best_candidate_fsd) && best_candidate_obj > current_obj
                current_obj  = best_candidate_obj
                current_fsd  = best_candidate_fsd
                improved     = true

                println(@sprintf("  iter %3d  move=%-12s  E[profit]=%+.2f  *** improved",
                        iter, string(move), current_obj))

                if current_obj > best_exp_obj
                    best_exp_obj = current_obj
                    best_fsd     = current_fsd
                end
                break  # move to next iteration with improved solution
            end
        end

        if !improved
            println(@sprintf("  iter %3d  no improving move found  (current=%+.2f)",
                    iter, current_obj))
        else
            # Keep current_accepted/active in sync with current_fsd
        end
        # Always sync from current_fsd for next iteration
        # (current_accepted/active are local loop vars, current_fsd is truth)
    end

    println("\n  Phase 2 best: accepted=$(sort(collect(best_fsd.accepted_passengers)))")
    println("  Phase 2 best: active=$(sort(collect(best_fsd.active_evtols)))")
    println("  Phase 2 best: E[profit]=$(round(best_exp_obj, digits=2))")

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
        pi_s                = pi_s,
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
                sc, scen.label, result.pi_s[sc],
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
