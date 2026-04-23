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
for file in sort(readdir(src_dir))
    endswith(file, ".jl") || continue
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
Maxtime = Int32(40) 

(best_obj, best_sol, iterations) = Heuristic(maxTurnaround, Maxtime, data, rt)

println("Heuristic ran $(iterations) iterations")
println("Best solution:")
println("Objective Value: $(best_obj)")
print_chromosome_table(best_sol)

assignments, scheduled = assign_passengersV2(best_sol, data, rt)

# println("Passenger Assignment")
print_assignments(assignments, data)

print_schedule_pretty(scheduled)
