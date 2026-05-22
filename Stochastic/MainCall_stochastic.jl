###############################################################################
# MainCall_stochastic.jl
#
# Entry point for the two-stage stochastic eVTOL heuristic.
# Place this file alongside the src/ folder and the inputData/ folder.
#
# Folder structure expected:
#   MainCall_stochastic.jl    ← this file
#   scenario_generation.jl    ← scenario parameters
#   src/
#       DataLoadFunc.jl
#       FeasibilityFunc.jl
#       InitialSolFunc.jl
#       PassAssignFunc.jl
#       FitnessFunc.jl
#       NeighborhoodFunc.jl
#       SANeighborhood.jl
#       HeuristicFunc.jl
#       StochasticHeuristic.jl  ← new file
#   inputData/
#       inputData.xlsx          (or inputDataMini.xlsx / inputDataGiant.xlsx)
#       Parameters.xlsx
###############################################################################

using XLSX
using DataFrames
using Random
using CSV
using Printf

src_dir = joinpath(@__DIR__, "..", "Heuristic", "src")

# Core heuristic files (same order as original MainCall.jl)
source_files = [
    "DataLoadFunc.jl",
    "FeasibilityFunc.jl",
    "InitialSolFunc.jl",
    "PassAssignFunc.jl",
    "FitnessFunc.jl",
    "NeighborhoodFunc.jl",
    "SANeighborhood.jl",
    "HeuristicFunc.jl"
]

for file in source_files
    include(joinpath(src_dir, file))
end

# Scenario generation lives next to this file (not inside src/)
include(joinpath(@__DIR__, "scenario_generation.jl"))
include(joinpath(@__DIR__, "StochasticHeuristic.jl"))

###############################################################################
# Configuration 
###############################################################################

excel_file     = joinpath(@__DIR__, "..", "inputData", "inputData.xlsx")
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")

maxTurnaround            = 100       # max turnaround time (minutes)
MaxTime                  = Int32(30) # first-stage heuristic time budget (s)
top_c                    = 4
sa_iters_per_scenario    = 30        # SA micro-improvement iters per scenario
                                     # increase for better quality, decrease for speed

###############################################################################
# Load data
###############################################################################

println("Loading data from: $excel_file")
data = load_data(excel_file, parameter_file)
println("  Vertiports : $(data.V)")
println("  eVTOLs     : $(data.N)")
println("  Passengers : $(data.A)")

println("Vertiports in V: $(data.V)")
println("Passenger origins/destinations: $([(data.op[a], data.dp[a]) for a in data.A])")

###############################################################################
# Generate scenarios
###############################################################################

time_per_km = derive_time_per_km(data)   # from StochasticHeuristic.jl
println("\nDerived time_per_km = $(round(time_per_km, digits=4)) min/km  " *
        "(v_air ≈ $(round(60/time_per_km, digits=1)) km/h)")

rt_s, e_s, S, pi_s = generate_scenarios(
    data.V, data.lat, data.lon, data.rt, data.e, time_per_km
)

###############################################################################
# Run two-stage stochastic heuristic
###############################################################################

result = stochastic_heuristic(
    data, rt_s, e_s, S, pi_s;
    maxTurnaround            = maxTurnaround,
    MaxTime                  = MaxTime,
    top_c                    = top_c,
    sa_iters_per_scenario    = sa_iters_per_scenario,
)

###############################################################################
# Print full results
###############################################################################

print_stochastic_result(result, data, rt_s, e_s)

###############################################################################
# Also print the deterministic best-solution details for comparison
###############################################################################

println("\n" * "="^72)
println("FIRST-STAGE SOLUTION — DETERMINISTIC DETAILS")
println("="^72)

Vmax = maximum(data.V)
rt_det = zeros(Vmax, Vmax)
for i in data.V, j in data.V
    rt_det[i, j] = data.rt[(i, j)]
end

feasible, battery_levels, _ = FeasibleBattery(
    result.first_stage_sol,
    Float32(data.bmax), Float32(data.bmid), Float32(data.bmin),
    data.dist, Float32(data.ec), Float32(data.battery_per_km),
    data.b_penalty
)
println("Battery levels (deterministic):")
for (p, levels) in enumerate(battery_levels)
    println("  Plane $p: $levels")
end

assignments_det, scheduled_det = assign_passengersV2(
    result.first_stage_sol, data, Int.(rt_det)
)
print_assignments(assignments_det, data)
print_schedule_pretty(scheduled_det)
