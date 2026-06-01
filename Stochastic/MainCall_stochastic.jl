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

<<<<<<< HEAD
MaxTime_1st          = Int32(5)
n_restarts           = 1
n_outer_iters        = 2
MaxTime_2nd_search   = Int32(3)
MaxTime_2nd_final    = Int32(10)

# ── Heuristic steering parameters ───────────────────────────────────────────
top_c        = 4
price_boost  = 100.0   # was 50 — needs to be dominant in route construction
hard_penalty = 100_000.0  # was 10,000 — needs to massively outweigh any other consideration
=======
# ── Phase 1: deterministic initialisation ───────────────────────────────────
MaxTime_1st = Int32(60)   # seconds per restart
n_restarts  = 3            # independent deterministic restarts

# ── Phase 2: outer stochastic search ────────────────────────────────────────
# Each iteration tries ALL candidates within each move type.
# With 13 scenarios × 10s each, one E[profit] evaluation ≈ 130s.
# Each iteration may evaluate O(|addable| + |accepted|) candidates.
# Keep n_outer_iters small (5-10) until runtime is acceptable.
n_outer_iters      = 10
MaxTime_2nd_search = Int32(20)  # seconds per scenario during search

# ── Phase 3: final reporting pass ───────────────────────────────────────────
MaxTime_2nd_final = Int32(60)   # longer budget for final routes

# ── Heuristic steering parameters ───────────────────────────────────────────
top_c        = 4
price_boost  = 50.0      # higher boost → heuristic more strongly prioritises
                          # accepted passenger OD pairs during route construction
hard_penalty = 10_000.0  # penalty per unserved accepted passenger
>>>>>>> 5dc14e8 (stochastic)

###############################################################################
# Load data and generate scenarios
###############################################################################

println("Loading data …")
<<<<<<< HEAD
excel_file     = joinpath(@__DIR__, "..", "inputData", "inputDataGiant.xlsx")
=======
excel_file     = joinpath(@__DIR__, "..", "inputData", "inputData.xlsx")
>>>>>>> 5dc14e8 (stochastic)
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
