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

excel_file     = joinpath(@__DIR__, "..", "inputData", "inputData.xlsx")
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")

maxTurnaround        = 100

# First-stage: how long to run HeuristicSA per restart to build tentative plan
MaxTime_1st          = Int32(30)

# Second-stage search: short budget per scenario during the search loop
# (trades solution quality for speed — keep 3-8s)
MaxTime_2nd_search   = Int32(5)

# Second-stage final: longer budget for the final reporting pass
MaxTime_2nd_final    = Int32(20)

top_c                = 4

# How much to boost accepted passenger fares in the second-stage heuristic
# so the route builder naturally visits their origins/destinations first.
# 50× means accepted passengers score 50× higher in k_BestRoutes.
price_boost          = 10.0

# Penalty per unserved accepted passenger in the second stage.
# Should be large enough to discourage leaving them unserved,
# but not so large it makes the objective uninformative.
# Rule of thumb: ~5-10× the typical fare for a single passenger.

# Number of independent first-stage restarts.
# Each restart builds a different tentative plan and extracts a different
# candidate (accepted_passengers, active_evtols) set.

hard_penalty         = 50_000.0   # well below 1_000_000 but strong enough
MaxTime_2nd_final    = Int32(30)  # more time to find feasible solutions
n_restarts           = 3         # more first-stage candidates

###############################################################################
# Load data and generate scenarios
###############################################################################

println("Loading data …")
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
)

print_stochastic_result(result, data)
