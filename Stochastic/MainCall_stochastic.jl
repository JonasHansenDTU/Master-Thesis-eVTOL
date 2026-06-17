###############################################################################
# MainCall_stochastic.jl
###############################################################################

using XLSX
using DataFrames
using Random
using CSV
using Printf

src_dir = joinpath(@__DIR__, "..", "Heuristic", "src")

source_files = [
    "DataLoadFunc.jl",
    "FeasibilityFunc.jl",
    "InitialSolFunc.jl",
    "PassAssignFunc.jl",
    "FitnessFunc.jl",
    "NeighborhoodFunc.jl",
    "SANeighborhood.jl",
    "HeuristicFunc.jl",
]

for file in source_files
    include(joinpath(src_dir, file))
end

include(joinpath(@__DIR__, "scenario_generation.jl"))
include(joinpath(@__DIR__, "StochasticHeuristic.jl"))

###############################################################################
# Configuration
###############################################################################

maxTurnaround = 100

MaxTime_1st        = Int32(90)   # MUST be enough to find the rich deterministic seed
n_restarts         = 2            # take the better of two seeds
n_outer_iters      = 8            # Phase 2 prunes/refines the seed (more exploration)
MaxTime_2nd_search = Int32(15)    # recourse quality during search (less noisy candidate scoring)
MaxTime_2nd_final  = Int32(30)    # final high-quality recourse pass
price_boost        = 8.0
hard_penalty       = 50_000.0
top_c              = 4            # top routes per base vertiport considered by constructor

###############################################################################
# Load data and generate scenarios
###############################################################################

println("Loading data …")
excel_file     = joinpath(@__DIR__, "..", "inputData", "inputDataGiant.xlsx")
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")

data = load_data(excel_file, parameter_file)
println("  Vertiports     : $(data.V)")
println("  eVTOLs         : $(length(data.N))")
println("  Passengers     : $(length(data.A))")
println("  Opening cost   : $(data.opening_cost)")
println("  p[a] (soft pen): $(data.p[first(data.A)])")
println("  cap_u (seats)  : $(data.cap_u)")

time_per_km = derive_time_per_km(data)
println("\nDerived time_per_km = $(round(time_per_km, digits=4)) min/km " *
        "(v_air ≈ $(round(60/time_per_km, digits=1)) km/h)")

rt_s, e_s, S, pi_s = generate_scenarios(
    data.V, data.lat, data.lon, data.rt, data.e, time_per_km
)

# ── Drop scenarios with zero probability in the active season ────────────────
# A scenario with pi = 0 contributes 0 to E[profit], so solving it is wasted
# work; removing it changes nothing in the objective. It also (correctly) means
# the robustness filter no longer requires passengers to be serveable in weather
# that cannot occur this season. rt_s / e_s are keyed by (sc,i,j) and need no
# pruning — only the scenario set S and probabilities pi_s drive the loops.
S_all          = collect(S)
S_active       = [sc for sc in S_all if pi_s[sc] > 0.0]
pi_s           = Dict(sc => pi_s[sc] for sc in S_active)
S              = S_active
dropped        = [sc for sc in S_all if !(sc in S_active)]

println("\nScenarios (active season): $(length(S)) of $(length(S_all))")
if !isempty(dropped)
    println("  Skipped zero-probability scenarios: $(dropped)")
end
println("Probabilities: $(round.([pi_s[sc] for sc in S], digits=4))")

###############################################################################
# Run two-stage stochastic heuristic
###############################################################################

###############################################################################
# Seed the global RNG so the first-stage solve is reproducible. Use the SAME
# value in Realised_scenario.jl so both produce the identical first-stage
# decision (the heuristic's SA draws from the global RNG).
###############################################################################
const FIRST_STAGE_SEED = 12345
Random.seed!(FIRST_STAGE_SEED)

result = stochastic_heuristic(
    data, rt_s, e_s, S, pi_s;
    maxTurnaround        = maxTurnaround,
    MaxTime_1st          = MaxTime_1st,
    MaxTime_2nd_search   = MaxTime_2nd_search,
    MaxTime_2nd_final    = MaxTime_2nd_final,
    top_c                = top_c,
    price_boost          = price_boost,
    hard_penalty         = hard_penalty,
    n_restarts           = n_restarts,
    n_outer_iters        = n_outer_iters,
)

print_stochastic_result(result, data)

###############################################################################
# Note: the stochastic-programming measures (VSS / EVPI) are computed by the
# separate driver RunMeasures.jl, so the plain two-stage run here stays fast and
# is unaffected by the (expensive) wait-and-see solves. Run RunMeasures.jl when
# you want EEV / WS / VSS / EVPI.
###############################################################################
