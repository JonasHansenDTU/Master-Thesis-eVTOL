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

    # active: set of eVTOL IDs (from data.N) that flew at least one leg.
    # enumerate(sol.planes) gives (position, plane) — position maps to data.N[position].
    active = Set{Int}(
        data.N[pos] for (pos, plane) in enumerate(sol.planes)
        if pos <= length(data.N) && plane.flightLegs > 0
    )

    # Robustness is evaluated against the FULL fleet, not just the currently
    # active eVTOLs. A passenger is robust if SOME eVTOL in the fleet could
    # serve it; the fleet-activation decision is handled separately in Phase 2.
    # Evaluating against only the active subset creates a chicken-and-egg trap
    # where few eVTOLs are active → few passengers pass → few eVTOLs needed.
    fleet_for_filter = Set{Int}(data.N)

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
                te_i      = Int(round(data.te))
                # Per-scenario operating window: severe weather relaxes ET, so
                # the filter is correspondingly less strict for those scenarios.
                _, ET_sc  = scenario_slack(data, sc)

                # Check 1: trip completes before this scenario's ET
                if data.dt[a] + trip_rt > ET_sc
                    robust = false; break
                end

                # Combined check: SOME SINGLE eVTOL n must do the whole chain in
                # this scenario — reach origin i from its base before dt[a], AND
                # after serving (arrive at j) plus turnaround te, return to one of
                # ITS OWN end vertiports before this scenario's ET. Evaluating
                # reach and return as one chain (rather than two independent
                # existence checks) prevents admitting a passenger that one eVTOL
                # can reach but only a different eVTOL could return from.
                # end_vp is keyed by vertiport ID, so use data.bv[n].
                can_serve = any(
                    (rt_mat_sc[data.bv[n], i] <= data.dt[a]) &&
                    (data.dt[a] + trip_rt + te_i +
                        minimum(rt_mat_sc[j, ev]
                                for ev in data.end_vp[data.bv[n]];
                                init = typemax(Int)) <= ET_sc)
                    for n in fleet_for_filter
                )
                if !can_serve
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

# Weather-dependent schedule slack. A scenario is "severe" if it has strong wind
# OR is very cold; severe scenarios get a single flat allowance (no stacking):
#   +10 min passenger waiting tolerance (w), +30 min operating window (ET).
# Returns (w_sc, ET_sc). The synthetic mean scenario (sc=0) and any out-of-range
# index keep baseline values.
function scenario_slack(data, sc::Int)
    w_sc  = Float64(data.w)
    ET_sc = Float64(data.ET)
    if sc >= 1 && sc <= length(SCENARIOS)
        scen        = SCENARIOS[sc]
        wind_speed  = sqrt(scen.wx^2 + scen.wy^2)
        strong_wind = wind_speed >= 45.0
        very_cold   = scen.phi >= 1.35
        if strong_wind || very_cold
            w_sc  += 10.0
            ET_sc += 30.0
        end
    end
    return w_sc, ET_sc
end

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

    # ─── Weather-dependent slack on waiting time (w) and operating window (ET) ──
    # See scenario_slack: severe scenarios (strong wind OR very cold) get a flat
    # +10 min on w and +30 min on ET; everything else keeps baseline.
    w_sc, ET_sc = scenario_slack(data, sc)

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
        w              = w_sc,
        ET             = ET_sc,
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

###############################################################################
# On-demand (Pool B) demand realizations
#
# Demand uncertainty: each potential on-demand group (data.B) independently
# "appears" in a given scenario with probability rho. The appearance is fixed
# per scenario with a deterministic seed so every run uses the SAME demand
# draws (reproducibility and comparability across runs/seasons).
#
# Each weather scenario carries ONE demand realization (the cheap, no-extra-
# scenarios design): scenario sc always has the same realized Pool B subset.
# This models weather and demand jointly through the existing scenario set.
#
# rho may be made weather-dependent by passing a function; by default it is a
# flat probability for every scenario.
###############################################################################

function make_demand_realizations(data, S; rho = 0.5, seed::Int = 20260101)
    realized = Dict{Int, Vector{Int}}()
    for sc in S
        # Per-scenario deterministic RNG so the realization is fixed and
        # reproducible, and independent across scenarios.
        rng = MersenneTwister(seed + sc)
        rho_sc = rho isa Function ? rho(sc) : rho
        appear = Int[]
        for b in data.B
            if rand(rng) < rho_sc
                push!(appear, b)
            end
        end
        realized[sc] = sort(appear)
    end
    return realized
end

function build_scenario_cache(data, rt_s, e_s, S; demand_realized = nothing)
    cache   = Dict{Int, NamedTuple}()
    rt_mats = Dict{Int, Matrix{Int}}()
    for sc in S
        sc_data = make_scenario_data(data, sc, rt_s, e_s)
        # Attach this scenario's realized on-demand (Pool B) groups so the
        # second stage knows which on-demand customers actually showed up.
        B_real = demand_realized === nothing ? Int[] : get(demand_realized, sc, Int[])
        sc_data = merge(sc_data, (B_realized = B_real,))
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

    # Set p[a] = hard_penalty for accepted passengers so the heuristic
    # is strongly penalised for missing them.
    # Set p[a] = original soft penalty for non-accepted passengers so the
    # heuristic still has incentive to serve them opportunistically.
    # This restores the original behaviour that was working before.
    new_p = copy(sc_data.p)
    for a in fsd.accepted_passengers
        new_p[a] = hard_penalty
    end

    # IMPORTANT: Construction_Heuristic2 uses eVTOL IDs directly as array indices
    # (planes[n] = ...) so N must always be the full fleet [1,2,3,4].
    # We cannot pass a subset like [3,4] because planes[3] would be out of bounds.
    # Non-committed eVTOLs are naturally discouraged because only committed eVTOLs
    # have their OD pairs fare-boosted, making them the most profitable to serve.

    return (
        infra          = sc_data.infra,
        pax            = sc_data.pax,
        plane          = sc_data.plane,
        V              = sc_data.V,
        A              = sc_data.A,
        N              = sc_data.N,  # must be full fleet — see comment above
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
        bmid           = sc_data.bmid,  # start at bmid (standing reserve charge)
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
        # Active-fleet restriction: the inner construction and SA neighborhoods
        # read this and build routes ONLY on these eVTOLs, so the second-stage
        # search optimises using exactly the committed fleet.
        active         = fsd.active_evtols,
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
                              hard_penalty::Float64 = 5_000.0,
                              active::Union{Nothing,Set{Int}} = nothing)

    # ── Active-fleet enforcement ────────────────────────────────────────────
    # Only first-stage activated eVTOLs may fly in the second stage (matching
    # the MIP, where routing of eVTOL n is gated by y[n]). Any eVTOL that flies
    # but is NOT in the active set makes this solution infeasible: we apply a
    # prohibitive penalty so the search avoids it, and treat its service as
    # invalid. (N must remain the full fleet for planes[n] indexing, so we
    # enforce activation here rather than by shrinking N.)
    fleet_violation = 0.0
    if active !== nothing
        for (idx, plane) in enumerate(sol.planes)
            if plane.flightLegs > 0 && !(idx in active)
                fleet_violation += 1_000_000.0
            end
        end
    end

    assignments, _ = assign_passengersV2(sol, sc_data, sc_rt_mat)
    served = Set{Int}(ass.group for ass in assignments)

    value = 0.0

    # Revenue from served passengers (true prices).
    # Must multiply by q[a] (group size) to match the MIP objective:
    #   q[a] * (fd[i,j]*(1-so[a]) + fs[i,j]*so[a])
    for ass in assignments
        a = ass.group
        i = sc_data.op[a]; j = sc_data.dp[a]
        value += sc_data.q[a] * (
                 sc_data.fd[(i,j)] * (1 - sc_data.so[a]) +
                 sc_data.fs[(i,j)] * sc_data.so[a])
    end

    # Operating cost for flown legs
    for plane in sol.planes
        for k in 1:plane.flightLegs
            value -= sc_data.c[(Int(plane.route[k]), Int(plane.route[k+1]))]
        end
    end

    # Physical feasibility penalties (battery, completion time, vertiport cap)
    # Start at bmid (standing reserve), cap at bmax. This matches the MIP's
    # initial battery assumption.
    P = FeasibilityCheck(
        Float32(sc_data.bmax), Float32(sc_data.bmid), Float32(sc_data.bmin),
        sc_data.dist, Float32(sc_data.ec), Float32(sc_data.battery_per_km),
        sol, sc_rt_mat,
        Int(round(sc_data.ET)), maximum(Int.(sc_data.T)),
        maximum(sc_data.V), sc_data.cap_v, sc_data.b_penalty
    )
    if sum(P) > 0
        println("    [feasibility] P = $P  " *
                "(battery=$(P[1]), completion_time=$(P[2]), vertiport_cap=$(P[3]))")
        if P[2] > 0
            # Show which eVTOL exceeds ET and by how much
            for (idx, plane) in enumerate(sol.planes)
                tt = 0
                for i in 1:plane.flightLegs
                    from = Int(plane.route[i]); to = Int(plane.route[i+1])
                    tt += Int(plane.turnaroundTime[i]) + sc_rt_mat[from, to]
                end
                if tt > Int(round(sc_data.ET))
                    println("      eVTOL pos $idx: route=$(Tuple(Int.(plane.route))) " *
                            "turnarounds=$(Int.(plane.turnaroundTime)) " *
                            "total_time=$tt > ET=$(Int(round(sc_data.ET)))")
                end
            end
        end
    end
    value -= sum(P)
    value -= fleet_violation   # prohibitive if a non-active eVTOL flew

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
                           accepted::Set{Int};
                           active::Union{Nothing,Set{Int}} = nothing)


    assignments, _ = assign_passengersV2(sol, sc_data, sc_rt_mat)
    served         = Set{Int}(ass.group for ass in assignments)
    still_unserved = Set{Int}()

    for a in accepted
        a in served && continue

        earliest = Int(round(sc_data.dt[a]))
        latest   = earliest + Int(round(sc_data.w))
        te       = Int(round(sc_data.te))
        fixed    = false

        # ── Attempt 1: nudge an EXISTING op→dp leg's turnaround into the window ──
        plane_indices = collect(enumerate(sol.planes))

        for (pidx, plane) in plane_indices
            plane.flightLegs == 0 && continue
            (active !== nothing && !(pidx in active)) && continue  # active fleet only
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

        # ── Attempt 2: serve via an IDLE active eVTOL based at the origin ───────
        # If no existing leg could be nudged, look for an activated but idle
        # eVTOL (flightLegs == 0) whose base is the passenger's origin vertiport.
        # Such an eVTOL can fly the single serving leg op→dp starting from the
        # ground, with a turnaround chosen so departure lands inside the time
        # window. This is safe because an idle eVTOL has no other legs to shift
        # and no other passengers to disrupt, so a single whole-solution
        # feasibility check is sufficient. (We deliberately do not insert into an
        # occupied route, which would shift later legs and could push other
        # passengers out of their windows or break ET / battery.)
        if !fixed
            for (pidx, plane) in plane_indices
                plane.flightLegs == 0 || continue            # must be idle
                (active !== nothing && !(pidx in active)) && continue  # active fleet only
                Int(plane.route[1]) == sc_data.op[a] || continue  # based at origin

                # Departure can be any time from 0; pick the earliest feasible
                # turnaround that lands inside [earliest, latest].
                new_ta  = clamp(earliest, 0, latest)
                new_dep = new_ta
                (earliest <= new_dep <= latest) || continue

                # Tentatively give the idle eVTOL the single serving leg.
                saved_route = copy(plane.route)
                saved_ta    = copy(plane.turnaroundTime)
                saved_legs  = plane.flightLegs

                plane.route          = Int32[Int32(sc_data.op[a]), Int32(sc_data.dp[a])]
                plane.turnaroundTime = Int32[Int32(new_ta)]
                plane.flightLegs     = Int32(1)

                P = FeasibilityCheck(
                    Float32(sc_data.bmax), Float32(sc_data.bmid), Float32(sc_data.bmin),
                    sc_data.dist, Float32(sc_data.ec), Float32(sc_data.battery_per_km),
                    sol, sc_rt_mat, Int(round(sc_data.ET)),
                    maximum(Int.(sc_data.T)), maximum(sc_data.V),
                    sc_data.cap_v, Float64(sc_data.b_penalty))

                feasible = all(p < 1_000_000 for p in P)

                if feasible
                    new_assign, _ = assign_passengersV2(sol, sc_data, sc_rt_mat)
                    if a in Set(ass.group for ass in new_assign)
                        fixed = true
                        break
                    end
                end

                # Revert if it didn't help or broke feasibility.
                plane.route          = saved_route
                plane.turnaroundTime = saved_ta
                plane.flightLegs     = saved_legs
            end
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

    # Second-stage candidate passengers =
    #   committed Pool A (must be served — carry the hard penalty)
    #   + this scenario's realized on-demand Pool B (opportunistic, no penalty).
    # Rejected Pool A passengers are intentionally EXCLUDED: they received a
    # "no" in the first stage and are not served in the second stage.
    B_real = hasproperty(sc_data, :B_realized) ? sc_data.B_realized : Int[]
    candidate_A = sort(collect(union(fsd.accepted_passengers, Set(B_real))))
    sc_data = merge(sc_data, (A = candidate_A,))

    ss_data = make_second_stage_data(sc_data, fsd;
                                      price_boost  = price_boost,
                                      hard_penalty = hard_penalty)

    # Run the heuristic. Because ss_data carries the active set, the construction
    # and SA neighbourhoods build routes ONLY on the committed eVTOLs, so the
    # search optimises directly over the activated fleet (no post-hoc stripping).
    _, sc_sol, _ = HeuristicSA(maxTurnaround, MaxTime_2nd, ss_data,
                                sc_rt_mat, top_c)

    active = fsd.active_evtols

    # Repair step: force-serve any accepted passengers still missed, using only
    # the committed (active) eVTOLs.
    repair_unserved!(sc_sol, sc_data, sc_rt_mat,
                     fsd.accepted_passengers; active = active)

    # ── Hard safety net ──────────────────────────────────────────────────────
    # The search is already confined to the active fleet, but as a guarantee we
    # strip any residual route on a non-active eVTOL before scoring. This ensures
    # the evaluated solution can NEVER use an un-activated aircraft, regardless of
    # any construction path. A Pool B group left unserved this way is fine (it is
    # opportunistic, no penalty); an accepted Pool A group would have been served
    # by repair on an active eVTOL above.
    for (idx, plane) in enumerate(sc_sol.planes)
        if !(idx in active) && plane.flightLegs > 0
            base = plane.route[1]
            plane.route          = Int32[Int32(base)]
            plane.turnaroundTime = Int32[]
            plane.flightLegs     = Int32(0)
        end
    end

    # Evaluate true profit on original scenario data, enforcing the active fleet.
    profit, assignments, accepted_served, n_unserved = second_stage_profit(
        sc_sol, sc_data, sc_rt_mat, fsd.accepted_passengers;
        hard_penalty = hard_penalty, active = active
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
        value += data.q[a] * (data.fd[(i,j)] * (1 - data.so[a]) + data.fs[(i,j)] * data.so[a])
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
    # Robustness is judged against the FULL fleet: a passenger is a candidate
    # if SOME eVTOL in the fleet could serve it in every scenario. The active
    # subset is not used here — fleet activation is a separate Phase 2 decision.
    # (The `active` argument is retained for call-site compatibility.)
    fleet = Set{Int}(data.N)
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
            te_i    = Int(round(data.te))
            _, ET_sc = scenario_slack(data, sc)   # severe weather relaxes ET
            if data.dt[a] + trip_rt > ET_sc
                robust = false; break
            end
            # Combined single-eVTOL chain: same eVTOL must reach the origin
            # before dt AND return to one of its own end vertiports before ET.
            can_serve = any(
                (rt_sc[data.bv[n], i] <= data.dt[a]) &&
                (data.dt[a] + trip_rt + te_i +
                    minimum(rt_sc[j, ev] for ev in data.end_vp[data.bv[n]];
                            init = typemax(Int)) <= ET_sc)
                for n in fleet
            )
            if !can_serve
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
                               n_outer_iters::Int        = 20,
                               demand_rho                = 0.5,
                               demand_seed::Int          = 20260101)

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
    demand_realized = make_demand_realizations(data, S;
                                               rho = demand_rho, seed = demand_seed)
    if !isempty(data.B)
        total_appear = sum(length(demand_realized[sc]) for sc in S)
        println("[Setup] On-demand pool B: $(length(data.B)) potential groups, " *
                "mean $(round(total_appear/length(S), digits=1)) appear per scenario " *
                "(rho = $(demand_rho)).")
    end
    scenario_cache, rt_mats = build_scenario_cache(data, rt_s, e_s, S;
                                                   demand_realized = demand_realized)
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

        # Cap how many candidates are tried per move type to control runtime.
        # With 13 scenarios × MaxTime_2nd_search per evaluation, each candidate
        # costs up to 13 × MaxTime_2nd_search seconds. At 5 candidates per move
        # type and 5 move types, one iteration costs at most:
        # 5 × 5 × 13 × MaxTime_2nd_search seconds.
        max_cand = 5

        for move in moves

            if move == :add
                for a in shuffle!(collect(addable))[1:min(max_cand, length(addable))]
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
                for a in shuffle!(collect(current_accepted))[1:min(max_cand, length(current_accepted))]
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
                for a_out in shuffle!(collect(current_accepted))[1:min(max_cand, length(current_accepted))]
                    for a_in in shuffle!(collect(addable))[1:min(2, length(addable))]
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

    # ─────────────────────────────────────────────────────────────────────────
    # Cleanup pass: drop any committed eVTOL that never actually flies.
    # The active set is seeded from the tentative plan and adjusted by FLEET±,
    # but the greedy search can leave a committed eVTOL that flies in NO scenario
    # (its opening cost is pure waste). We evaluate the current best, find which
    # active eVTOLs carry no leg in any scenario, and drop them — this can only
    # save opening cost without losing service, so it never worsens E[profit].
    # ─────────────────────────────────────────────────────────────────────────
    println("\n[Cleanup] Checking for committed eVTOLs that never fly …")
    _, cleanup_results = expected_objective(
        best_fsd, scenario_cache, rt_mats, S, pi_s;
        maxTurnaround = maxTurnaround, MaxTime_2nd = MaxTime_2nd_search,
        top_c = top_c, price_boost = price_boost,
        hard_penalty = hard_penalty, full_output = true,
    )
    flew = Set{Int}()
    for (sc, r) in cleanup_results
        for (idx, plane) in enumerate(r.sol.planes)
            plane.flightLegs > 0 && push!(flew, idx)
        end
    end
    never_fly = setdiff(best_fsd.active_evtols, flew)
    if !isempty(never_fly)
        trimmed_active = setdiff(best_fsd.active_evtols, never_fly)
        if !isempty(trimmed_active)
            trimmed_fsd = make_fsd(best_fsd.accepted_passengers,
                                   Set{Int}(trimmed_active), data, best_tent_sol)
            trimmed_obj, _ = expected_objective(
                trimmed_fsd, scenario_cache, rt_mats, S, pi_s;
                maxTurnaround = maxTurnaround, MaxTime_2nd = MaxTime_2nd_search,
                top_c = top_c, price_boost = price_boost,
                hard_penalty = hard_penalty, full_output = false,
            )
            println("  Dropping never-flying eVTOLs $(sort(collect(never_fly))) " *
                    "(saves opening cost). E[profit] $(round(best_exp_obj,digits=2)) " *
                    "→ $(round(trimmed_obj,digits=2))")
            best_fsd     = trimmed_fsd
            best_exp_obj = trimmed_obj
        end
    else
        println("  None — every committed eVTOL flies in at least one scenario.")
    end

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

    # ─────────────────────────────────────────────────────────────────────────
    # On-demand (Pool B) demand realization per scenario: which on-demand groups
    # ACTUALLY APPEARED in each scenario (the demand-uncertainty realization),
    # shown alongside the scenario probability. "appeared" is independent of
    # whether they were ultimately served — it is the realized demand the second
    # stage had the option to serve. Pool B IDs shown with the +1000 offset removed.
    # ─────────────────────────────────────────────────────────────────────────
    if !isempty(data.B)
        println("\n" * "="^72)
        println("ON-DEMAND DEMAND REALIZATION (Pool B groups present per scenario)")
        println("="^72)
        println(@sprintf("  %-4s  %-18s  %8s  %7s  %s",
                "sc", "label", "prob", "#appear", "on-demand groups that appeared"))
        println("  " * "-"^78)
        for sc in sort(collect(keys(result.scenario_results)))
            scen   = SCENARIOS[sc]
            r      = result.scenario_results[sc]
            b_real = hasproperty(r.sc_data, :B_realized) ? r.sc_data.B_realized : Int[]
            b_str  = isempty(b_real) ? "—" :
                     join([g > 1000 ? g - 1000 : g for g in sort(b_real)], ", ")
            println(@sprintf("  %-4d  %-18s  %8.4f  %7d  [%s]",
                    sc, scen.label, result.pi_s[sc], length(b_real), b_str))
        end
        println("  " * "-"^78)
        println("  (Pool B IDs shown with the +1000 offset removed; appearance is the")
        println("   realized demand — see the service table for which were actually served.)")
    end

    # ─────────────────────────────────────────────────────────────────────────
    # Weather slack actually in effect per scenario: shows baseline vs relaxed
    # waiting time (w) and operating window (ET). Severe scenarios (strong wind
    # OR very cold) receive +10 w and +30 ET.
    # ─────────────────────────────────────────────────────────────────────────
    println("\n" * "="^72)
    println("WEATHER SLACK: allowed vs actually used")
    println("="^72)
    println(@sprintf("  %-4s  %-18s  %4s %4s  %9s  %4s %4s  %s",
            "sc", "label", "w", "ET", "relaxed?", "useW", "useET", "binding?"))
    println("  " * "-"^74)
    base_w  = Float64(data.w)
    base_ET = Float64(data.ET)
    for sc in sort(collect(keys(result.scenario_results)))
        scen        = SCENARIOS[sc]
        w_sc, ET_sc = scenario_slack(data, sc)
        relaxed     = (w_sc > base_w) || (ET_sc > base_ET)
        r           = result.scenario_results[sc]

        # Rebuild the realized schedule for this scenario to measure what was
        # actually used: the latest arrival time (realized completion) and the
        # largest passenger wait (boarding dep minus earliest-ready dt).
        Vmax = maximum(r.sc_data.V)
        rt_mat_sc = zeros(Int, Vmax, Vmax)
        for (ij, v) in r.sc_data.rt
            rt_mat_sc[ij[1], ij[2]] = v
        end
        sched = build_scheduled_legs(r.sol, rt_mat_sc, Int(round(r.sc_data.cap_u)))
        used_ET = isempty(sched) ? 0 : maximum(l.arr for l in sched)

        used_w = 0
        for ass in r.assignments
            isempty(ass.legs) && continue
            board_leg_idx = first(ass.legs)
            board = findfirst(l -> l.plane == ass.plane && l.leg_index == board_leg_idx, sched)
            board === nothing && continue
            wait = sched[board].dep - Int(round(r.sc_data.dt[ass.group]))
            used_w = max(used_w, wait)
        end

        # "binding" if the realized usage actually exceeded the baseline limits,
        # i.e. the relaxation was not just available but actually needed.
        binds_ET = used_ET > base_ET
        binds_w  = used_w  > base_w
        bind_str = (binds_ET || binds_w) ?
                   (binds_ET && binds_w ? "ET & w" : (binds_ET ? "ET" : "w")) : "—"

        println(@sprintf("  %-4d  %-18s  %4.0f %4.0f  %9s  %4d %4d  %s",
                sc, scen.label, w_sc, ET_sc,
                relaxed ? "yes" : "—", used_w, used_ET, bind_str))
    end
    println("  " * "-"^74)
    println(@sprintf("  Baseline w = %.0f, ET = %.0f   |   Severe adds +10 w, +30 ET",
            base_w, base_ET))
    println("  useW / useET = largest wait and latest arrival actually realized.")
    println("  binding = realized usage exceeded the baseline limit (slack was needed).")

    # ─────────────────────────────────────────────────────────────────────────
    # Per-scenario passenger service breakdown, split by pool:
    #   commitA  = first-stage committed (Pool A) groups served this scenario
    #   missA    = committed (Pool A) groups NOT served (incur the penalty)
    #   onDemB   = on-demand (Pool B) groups served opportunistically
    # Pool B groups are identified by membership in data.B.
    # ─────────────────────────────────────────────────────────────────────────
    poolB = Set{Int}(data.B)
    println("\n" * "="^72)
    println("SECOND-STAGE PASSENGER SERVICE (committed Pool A  vs  on-demand Pool B)")
    println("="^72)
    println(@sprintf("  %-4s  %-18s  %7s  %6s  %7s  %s",
            "sc", "label", "commitA", "missA", "onDemB", "on-demand groups served"))
    println("  " * "-"^80)
    accepted = result.accepted_passengers
    total_B = 0
    total_missA = 0
    for sc in sort(collect(keys(result.scenario_results)))
        scen = SCENARIOS[sc]
        r    = result.scenario_results[sc]
        served = Set(ass.group for ass in r.assignments)

        commitA_served = sort([a for a in accepted if a in served])
        commitA_missed = sort([a for a in accepted if !(a in served)])
        onDemB_served  = sort([g for g in served if g in poolB])

        total_B     += length(onDemB_served)
        total_missA += length(commitA_missed)

        # Show Pool B IDs with the 1000-offset stripped for readability.
        b_str = isempty(onDemB_served) ? "—" :
                join([g > 1000 ? g - 1000 : g for g in onDemB_served], ", ")

        println(@sprintf("  %-4d  %-18s  %7d  %6d  %7d  [%s]",
                sc, scen.label,
                length(commitA_served), length(commitA_missed),
                length(onDemB_served), b_str))
    end
    nsc = length(result.scenario_results)
    println("  " * "-"^80)
    println(@sprintf("  Mean on-demand (Pool B) groups served per scenario : %.2f", total_B / nsc))
    println(@sprintf("  Total committed (Pool A) misses across scenarios   : %d", total_missA))
    println("  (Pool B IDs shown with the +1000 offset removed.)")

    # ─────────────────────────────────────────────────────────────────────────
    # Key performance indicators
    # ─────────────────────────────────────────────────────────────────────────
    profits = [result.scenario_results[sc].profit
               for sc in sort(collect(keys(result.scenario_results)))]
    probs   = [result.pi_s[sc]
               for sc in sort(collect(keys(result.scenario_results)))]
    exp_ss  = sum(p*v for (p,v) in zip(probs, profits))
    # probability-weighted variance / std of second-stage profit
    var_ss  = sum(p*(v - exp_ss)^2 for (p,v) in zip(probs, profits))
    std_ss  = sqrt(max(0.0, var_ss))
    min_profit = minimum(profits)
    max_profit = maximum(profits)
    n_negative = count(<(0), profits)

    # mean total groups served per scenario (Pool A committed-served + Pool B)
    mean_served = sum(probs[k] * length(Set(ass.group for ass in
                      result.scenario_results[sc].assignments))
                      for (k, sc) in enumerate(sort(collect(keys(result.scenario_results)))))

    # mean on-demand (Pool B) groups served per scenario, probability-weighted
    mean_onDemB = sum(probs[k] * count(g -> g in poolB,
                      Set(ass.group for ass in result.scenario_results[sc].assignments))
                      for (k, sc) in enumerate(sort(collect(keys(result.scenario_results)))))

    println("\n" * "="^72)
    println("KEY PERFORMANCE INDICATORS")
    println("="^72)
    @printf("  E[total profit]            (RP)        : %+12.2f\n", result.expected_obj)
    @printf("  E[second-stage profit]                 : %+12.2f\n", exp_ss)
    @printf("  First-stage (opening) cost             : %12.2f\n", result.first_stage_cost)
    @printf("  Std-dev of 2nd-stage profit            : %12.2f\n", std_ss)
    @printf("  Min / Max scenario profit              : %+10.2f / %+.2f\n", min_profit, max_profit)
    @printf("  Scenarios with negative profit         : %12d / %d\n", n_negative, length(profits))
    @printf("  Committed (Pool A) passengers          : %12d\n", length(result.accepted_passengers))
    @printf("  On-demand pool (Pool B) size           : %12d\n", length(poolB))
    @printf("  Active eVTOLs                          : %12d\n", length(result.active_evtols))
    @printf("  Mean groups served / scenario (total)  : %12.2f\n", mean_served)
    @printf("  Mean on-demand (Pool B) served / scen. : %12.2f\n", mean_onDemB)
    println("="^72)

    println("\n" * "="^72)
    println("FINAL SUMMARY")
    println("="^72)
    @printf("  First-stage cost (opening)     = %+.2f\n", -result.first_stage_cost)
    @printf("  E[second-stage profit]         = %+.2f\n", exp_ss)
    @printf("  E[total profit]  (E[profit])   = %+.2f\n", result.expected_obj)
    println("="^72)

    # ─────────────────────────────────────────────────────────────────────────
    # Per-scenario second-stage routing audit. Prints the ACTUAL routes flown in
    # each scenario (the binding solution, not the tentative seed), and for each
    # leg the passengers carried with a seat-capacity check. Flags any eVTOL that
    # flies but is NOT in the committed active set (active-set leak), and any leg
    # whose carried passengers exceed cap_u (capacity-accounting bug).
    # ─────────────────────────────────────────────────────────────────────────
    print_scenario_routing_audit(result, data)
end

function print_scenario_routing_audit(result, data)
    active = Set(result.active_evtols)
    cap_u  = Int(round(data.cap_u))
    println("\n" * "="^72)
    println("PER-SCENARIO SECOND-STAGE ROUTING AUDIT")
    println("="^72)
    println("  Active (committed) eVTOLs: $(sort(collect(active)))")
    println("  Seat capacity per leg (cap_u): $cap_u")

    for sc in sort(collect(keys(result.scenario_results)))
        scen = SCENARIOS[sc]
        r    = result.scenario_results[sc]
        println("\n  " * "-"^68)
        println("  Scenario $sc — $(scen.label)")
        println("  " * "-"^68)
        println("  Routes actually flown (second-stage solution):")
        print_chromosome_table(r.sol)
        println("  Passenger-carrying legs and capacity check:")

        # Group assignments by plane and leg index.
        # ass.legs lists the operation indices on which group ass.group is carried.
        leg_pax = Dict{Tuple{Int,Int}, Vector{Int}}()   # (plane, leg) -> [groups]
        planes_used = Set{Int}()
        for ass in r.assignments
            push!(planes_used, ass.plane)
            for lg in ass.legs
                push!(get!(leg_pax, (ass.plane, lg), Int[]), ass.group)
            end
        end

        # Report each plane that actually flies in this scenario.
        for n in sort(collect(planes_used))
            leak = n in active ? "" : "  ⚠ NOT IN ACTIVE SET (leak?)"
            println("    eVTOL $n$leak")
            legs_here = sort([lg for (pl,lg) in keys(leg_pax) if pl == n])
            if isempty(legs_here)
                println("      (flies but carries no committed/served group)")
            end
            for lg in legs_here
                groups = sort(leg_pax[(n,lg)])
                used   = sum(data.q[g] for g in groups)
                over   = used > cap_u ? "  ⚠ EXCEEDS cap_u=$cap_u" : ""
                gstr   = join([g > 1000 ? "B$(g-1000)" : "A$g" for g in groups], ", ")
                println(@sprintf("      leg %-2d : seats %d/%d  [%s]%s",
                        lg, used, cap_u, gstr, over))
            end
        end
        if isempty(planes_used)
            println("    (no eVTOL carries any group in this scenario)")
        end
    end
    println("\n" * "="^72)
    println("  Legend: A# = committed Pool A group, B# = on-demand Pool B group.")
    println("  ⚠ flags an eVTOL flying outside the active set, or a leg over capacity.")
    println("="^72)
end
