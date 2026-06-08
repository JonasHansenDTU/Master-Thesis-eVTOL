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
n_outer_iters      = 4            # Phase 2 prunes/refines the seed
MaxTime_2nd_search = Int32(8)     # recourse quality during search
MaxTime_2nd_final  = Int32(30)    # final high-quality recourse pass
price_boost        = 10.0
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

println("\nScenarios: $(length(S))")
println("Probabilities: $(round.(values(pi_s) |> collect, digits=4))")

###############################################################################
# Run two-stage stochastic heuristic
###############################################################################

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
# Build the scenario cache (needed by the VSS/EVPI analysis). This mirrors what
# stochastic_heuristic builds internally.
###############################################################################
scenario_cache, rt_mats = build_scenario_cache(data, rt_s, e_s, S)

###############################################################################
# Stochastic-programming measures: VSS and EVPI
#
# RP is the stochastic solution's E[profit] (result.expected_obj).
# EEV solves the mean scenario, fixes that decision, evaluates across scenarios.
# WS solves each scenario with perfect foresight and weights by probability.
# Uses a generous deterministic budget so WS/EEV are well approximated.
###############################################################################
# include(joinpath(@__DIR__, "compute_vss_evpi.jl"))

# RP = result.expected_obj

# run_vss_evpi(
#     data, rt_s, e_s, S, pi_s, scenario_cache, rt_mats, RP;
#     maxTurnaround = maxTurnaround,
#     MaxTime_det   = Int32(90),     # generous per-scenario deterministic budget
#     MaxTime_2nd   = MaxTime_2nd_final,
#     top_c         = top_c,
#     price_boost   = price_boost,
#     hard_penalty  = hard_penalty,
# )
