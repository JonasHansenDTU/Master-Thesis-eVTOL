using XLSX
using DataFrames
using Random
using JuMP
using Gurobi
using CSV
using MathOptInterface
using Printf
const MOI = MathOptInterface


src_dir = joinpath(@__DIR__, "src")
source_files = [
    "DataLoadFunc.jl",
    "FeasibilityFunc.jl",
    "InitialSolFunc.jl",
    "PassAssignFunc.jl",
    "FitnessFunc.jl",
    "NeighborhoodFunc.jl",
    "HeuristicFunc.jl",
]

for file in source_files
    include(joinpath(src_dir, file))
end




excel_file = joinpath(@__DIR__, "inputData.xlsx")
data = load_data(excel_file)

Vmax = maximum(data.V)
rt = zeros(Int, Vmax, Vmax)
for i in data.V, j in data.V
    rt[i, j] = data.rt[(i, j)]
end


maxTurnaround = 100
Maxtime = Int32(30) 

(best_obj, best_sol, iterations) = Heuristic(maxTurnaround, Maxtime, data, rt)

println("Heuristic ran $(iterations) iterations")
println("Best solution:")
println("Objective Value: $(best_obj)")
print_chromosome_table(best_sol)

# Compute and print battery levels
feasible, battery_levels = FeasibleBattery(best_sol, Float32(data.bmax), Float32(data.bmin), data.dist, Float32(data.ec), Float32(data.battery_per_km))
println("Battery Levels:")
for (p, levels) in enumerate(battery_levels)
    println("Plane $p: $(levels)")
end

assignments, scheduled = assign_passengersV2(best_sol, data, rt)

# println("Passenger Assignment")
print_assignments(assignments, data)

print_schedule_pretty(scheduled)

# Export solution snapshots for visualization
snapshots = export_solution_snapshots(best_sol, scheduled, assignments, battery_levels, data)
