###############################################################################
# Realised_scenario.jl
#
# OUT-OF-SAMPLE EVALUATION of the two-stage solution.
#
# The first stage commits Pool A passengers and activates eVTOLs by hedging
# across the 13 planning scenarios. This file tests how that fixed commitment
# performs on "real" operating days that are drawn independently and need NOT be
# one of the 13 planning scenarios:
#
#   1. Run the first stage once  → fix (accepted passengers, active eVTOLs).
#   2. For each of N_DAYS:
#        a. draw a realised weather (wind direction + temperature) from the
#           season distribution — any wind×temp combination, not just the 13;
#        b. draw a realised on-demand (Pool B) set (same ρ as planning, fresh);
#        c. build that day's data (weather-adjusted travel/energy, slack);
#        d. solve the second stage ONCE with the commitment fixed;
#        e. record realised profit and service.
#   3. Report the distribution of realised profit and service over the N_DAYS,
#      and compare its mean against the in-sample RP.
#
# Set N_DAYS = 1 for a single illustrative day, or e.g. 100 for a proper
# out-of-sample distribution.
###############################################################################

using XLSX, DataFrames, Random, CSV, Printf, Statistics

###############################################################################
# CONFIGURATION
###############################################################################

const ACTIVE_SEASON = :Sommer      # :Foraar, :Sommer, :Efteraar, :Vinter
const N_DAYS        = 100          # number of realised days to simulate
const REALISED_SEED = 4242         # seed for the realised-day draws (reproducible)

# Temperature penalty (phi) levels and seasonal mix (DMI-DATA.xlsx Sheet2)
const PHI_MAP = [1.00, 1.20, 1.35]
const TEMP_PROBS = Dict(
    :Foraar   => (0.825, 0.153, 0.022),
    :Sommer   => (1.000, 0.000, 0.000),
    :Efteraar => (0.934, 0.055, 0.011),
    :Vinter   => (0.629, 0.293, 0.078),
)

# Wind-direction options (wx, wy in km/h) and shares (DMI TR12-19 wind roses)
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

# Two-stage / second-stage parameters — keep identical to MainCall_stochastic.jl
maxTurnaround      = 100
MaxTime_1st        = Int32(90)
n_restarts         = 2
n_outer_iters      = 4
MaxTime_2nd_search = Int32(8)
MaxTime_2nd_final  = Int32(30)
price_boost        = 8.0
hard_penalty       = 50_000.0
top_c              = 4
demand_rho         = 0.5
demand_seed        = 20260101

###############################################################################
# Suppress the verbose internal logging of HeuristicSA / the pool dump, so only
# this driver's own summary is printed. Runs the given function with stdout
# redirected to a null sink, then restores stdout and returns the result.
###############################################################################
function quiet(f)
    # Redirect to devnull (NOT a pipe): the heuristic prints a large volume of
    # logs, and redirecting to an unread pipe deadlocks once its buffer fills.
    # devnull is already an open IO sink, so redirect straight to it.
    redirect_stdout(f, devnull)
end

###############################################################################
# Load sources
###############################################################################

src_dir = joinpath(@__DIR__, "..", "Heuristic", "src")
for f in ["DataLoadFunc.jl","FeasibilityFunc.jl","InitialSolFunc.jl",
          "PassAssignFunc.jl","FitnessFunc.jl","NeighborhoodFunc.jl",
          "SANeighborhood.jl","HeuristicFunc.jl"]
    include(joinpath(src_dir, f))
end
include(joinpath(@__DIR__, "scenario_generation.jl"))
include(joinpath(@__DIR__, "StochasticHeuristic.jl"))

###############################################################################
# Sampling: one realised weather day
###############################################################################

# Robust categorical sampler: returns the index whose cumulative prob first
# exceeds u; guards against floating-point overshoot.
function sample_index(weights::Vector{Float64}, u::Float64)
    c = 0.0
    for (i, w) in enumerate(weights)
        c += w
        if u <= c
            return i
        end
    end
    return length(weights)
end

function draw_weather(season::Symbol, rng)
    tp     = collect(TEMP_PROBS[season])
    t_lvl  = sample_index(tp ./ sum(tp), rand(rng))
    phi    = PHI_MAP[t_lvl]
    t_lab  = ["AboveZero","Cold","VeryCold"][t_lvl]

    wsh    = [w.share for w in WIND_OPTIONS]
    w_idx  = sample_index(wsh ./ sum(wsh), rand(rng))
    w      = WIND_OPTIONS[w_idx]

    return (wx=w.wx, wy=w.wy, phi=phi, label="$(w.key) / $(t_lab)")
end

function draw_demand(data, rng; rho = 0.5)
    appear = Int[]
    for b in data.B
        if rand(rng) < rho
            push!(appear, b)
        end
    end
    return sort(appear)
end

###############################################################################
# Build one realised day's sc_data from a drawn weather realisation.
# Mirrors the physics in scenario_generation.jl and applies the same
# weather-dependent w / ET slack (severe = strong wind OR very cold).
###############################################################################

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

    # Weather-dependent slack: severe = strong wind (speed >= 45) OR very cold.
    wind_speed  = sqrt(weather.wx^2 + weather.wy^2)
    severe      = (wind_speed >= 45.0) || (weather.phi >= 1.35)
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
# Load data, scenarios, and run the FIRST STAGE once
###############################################################################

println("Loading data …")
excel_file     = joinpath(@__DIR__, "..", "inputData", "inputDataGiant.xlsx")
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")
data        = load_data(excel_file, parameter_file)
time_per_km = derive_time_per_km(data)

rt_s, e_s, S, pi_s = generate_scenarios(data.V, data.lat, data.lon, data.rt, data.e, time_per_km)
S_active = [sc for sc in collect(S) if pi_s[sc] > 0.0]
pi_s = Dict(sc => pi_s[sc] for sc in S_active); S = S_active

###############################################################################
# Seed the global RNG identically to MainCall_stochastic.jl so the first-stage
# decision here is the SAME one MainCall produces (the heuristic's SA draws from
# the global RNG). This is seeded BEFORE the first-stage solve; the realised-day
# draws use their own separate RNG (REALISED_SEED) afterwards.
###############################################################################
const FIRST_STAGE_SEED = 12345
Random.seed!(FIRST_STAGE_SEED)

println("\nRunning the two-stage heuristic to fix the first-stage decision …")
result = quiet() do
    stochastic_heuristic(
        data, rt_s, e_s, S, pi_s;
        maxTurnaround = maxTurnaround, MaxTime_1st = MaxTime_1st,
        MaxTime_2nd_search = MaxTime_2nd_search, MaxTime_2nd_final = MaxTime_2nd_final,
        top_c = top_c, price_boost = price_boost, hard_penalty = hard_penalty,
        n_restarts = n_restarts, n_outer_iters = n_outer_iters,
        demand_rho = demand_rho, demand_seed = demand_seed,
    )
end
RP  = result.expected_obj
fsd = FirstStageDecision(result.accepted_passengers, result.active_evtols,
                         result.tentative_sol, result.first_stage_cost)

println("\nFirst-stage decision fixed:")
println("  Committed Pool A : $(sort(collect(result.accepted_passengers)))")
println("  Active eVTOLs    : $(sort(collect(result.active_evtols)))")
println("  In-sample RP     : $(round(RP, digits=2))")

###############################################################################
# Simulate N_DAYS realised days with the commitment fixed
###############################################################################

rng = MersenneTwister(REALISED_SEED)
day_profit   = Float64[]
day_missA    = Int[]
day_onDemB   = Int[]
day_severe   = Int[]
day_records  = Vector{NamedTuple}(undef, 0)
day_full     = Any[]                       # full second-stage result per day
day_Breal    = Vector{Vector{Int}}(undef, 0)  # realized Pool B per day

println("\nSimulating $(N_DAYS) realised day(s) …")
for d in 1:N_DAYS
    weather = draw_weather(ACTIVE_SEASON, rng)
    B_real  = draw_demand(data, rng; rho = demand_rho)
    sc_data, rt_mat, severe = build_realised_day_data(data, weather, B_real, time_per_km)

    r = quiet() do
        solve_second_stage(fsd, sc_data, rt_mat;
                           maxTurnaround = maxTurnaround,
                           MaxTime_2nd   = MaxTime_2nd_final,
                           top_c = top_c, price_boost = price_boost,
                           hard_penalty = hard_penalty)
    end

    served  = Set(ass.group for ass in r.assignments)
    poolB   = Set(data.B)
    missA   = count(a -> !(a in served), collect(result.accepted_passengers))
    onDemB  = count(g -> g in poolB, served)

    push!(day_profit, r.profit)
    push!(day_missA, missA)
    push!(day_onDemB, onDemB)
    push!(day_severe, severe ? 1 : 0)
    push!(day_full, r)
    push!(day_Breal, B_real)
    push!(day_records, (day=d, label=weather.label, profit=r.profit,
                        missA=missA, onDemB=onDemB, appearedB=length(B_real)))
end

###############################################################################
# Report
###############################################################################

if N_DAYS == 1
    rec    = day_records[1]
    r1     = day_full[1]                 # the full second-stage result for the day
    poolB  = Set(data.B)
    served = Set(ass.group for ass in r1.assignments)
    B_real = day_Breal[1]

    committed     = sort(collect(result.accepted_passengers))
    commit_served = [a for a in committed if a in served]
    commit_missed = [a for a in committed if !(a in served)]
    onDem_served  = sort([g for g in served if g in poolB])

    strip_b(v) = [g > 1000 ? g - 1000 : g for g in v]

    println("\n" * "="^60)
    println("REALISED DAY (out-of-sample)")
    println("="^60)
    println("  Weather            : $(rec.label)")
    println("  Severe slack used  : $(day_severe[1] == 1 ? "yes (w+10, ET+30)" : "no")")
    println("  Realised profit    : $(round(rec.profit, digits=2))")
    println("  In-sample RP       : $(round(RP, digits=2))")
    println("-"^60)
    println("  PASSENGER SERVICE")
    println("    Committed (Pool A) served  : $(commit_served)")
    println("    Committed (Pool A) MISSED  : $(isempty(commit_missed) ? "none" : commit_missed)")
    println("    On-demand (Pool B) appeared: $(strip_b(sort(B_real)))")
    println("    On-demand (Pool B) served  : $(isempty(onDem_served) ? "none" : strip_b(onDem_served))")
    println("    On-demand capture          : $(length(onDem_served)) of $(length(B_real)) appeared")
    println("-"^60)
    println("  ROUTES FLOWN THIS DAY")
    print_chromosome_table(r1.sol)
    println("="^60)
else
    println("\n" * "="^72)
    println("OUT-OF-SAMPLE EVALUATION over $(N_DAYS) realised days  (season: $(ACTIVE_SEASON))")
    println("="^72)
    @printf("  Mean realised profit        : %+12.2f\n", mean(day_profit))
    @printf("  Std-dev realised profit     : %12.2f\n", std(day_profit))
    @printf("  Min / Max realised profit   : %+10.2f / %+.2f\n",
            minimum(day_profit), maximum(day_profit))
    @printf("  In-sample RP (for compare)  : %+12.2f\n", RP)
    @printf("  Days with a committed miss  : %12d / %d\n", count(>(0), day_missA), N_DAYS)
    @printf("  Mean committed misses/day   : %12.3f\n", mean(day_missA))
    @printf("  Mean on-demand served/day   : %12.2f\n", mean(day_onDemB))
    # On-demand capture rate = served / appeared, averaged over days that had
    # any on-demand appear (days with zero appearances are excluded to avoid 0/0).
    appeared = [length(b) for b in day_Breal]
    capture_days = [day_onDemB[i] / appeared[i] for i in 1:N_DAYS if appeared[i] > 0]
    if !isempty(capture_days)
        @printf("  Mean on-demand capture rate : %11.1f%%  (served / appeared)\n",
                100 * mean(capture_days))
    end
    @printf("  Severe-weather days drawn   : %12d / %d\n", sum(day_severe), N_DAYS)
    println("="^72)
    # Optional: save per-day records to CSV for plotting in the thesis.
    out = DataFrame(day_records)
    CSV.write(joinpath(@__DIR__, "realised_days.csv"), out)
    println("  Per-day records written to realised_days.csv")
end
