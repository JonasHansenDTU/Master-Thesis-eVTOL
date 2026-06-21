using XLSX
using DataFrames
using Random
using JuMP
using Gurobi
using CSV
using Statistics
using MathOptInterface
using Printf
const MOI = MathOptInterface

src_dir = joinpath(@__DIR__,"..",  "src")
source_files = [
    "DataLoadFunc.jl",
    "FeasibilityFunc.jl",
    "InitialSolFunc.jl",
    "PassAssignFunc.jl",
    "KPIFunc.jl",
    "FitnessFunc.jl",
    "NeighborhoodFunc.jl",
    "SANeighborhood.jl",
    "HeuristicFunc.jl",
]

for file in source_files
    include(joinpath(src_dir, file))
end

excel_file     = joinpath(@__DIR__, "..", "..", "inputData", "Experiments/inputDataEx2_1_2.xlsx")
parameter_file = joinpath(@__DIR__, "..", "..", "inputData", "Parameters.xlsx")
data = load_data(excel_file, parameter_file)

Vmax = maximum(data.V)
rt = zeros(Vmax, Vmax)
rt = zeros(Vmax, Vmax)
for i in data.V, j in data.V
    rt[i,j] = data.rt[(i,j)]
    rt[i,j] = data.rt[(i,j)]
end

maxTurnaround = Int64(data.ET)
Maxtime = Int32(10)
top_c = 4
Optimum = 4923.38
# Random.seed!(1234)

LoopIter = 5

obj_vals = zeros(LoopIter)
iter_vals = zeros(LoopIter)
ttg_vals = zeros(LoopIter)
tto_vals = zeros(LoopIter)
ttb_vals = zeros(LoopIter)
profit_vals = zeros(LoopIter)

for i in 1:LoopIter
    (obj_vals[i], best_sols, iter_vals[i], ttg_vals[i], tto_vals[i], ttb_vals[i], profit_vals[i]) =
        HeuristicSA(maxTurnaround, Maxtime, data, rt, top_c; optimal_obj = Float32(Optimum))
end


avg_obj = mean(obj_vals)
avg_iter = mean(iter_vals)
if Inf in ttg_vals
    avg_ttg = Inf32
else 
    avg_ttg = mean(ttg_vals)
end
if Inf in tto_vals
    avg_tto = Inf32
else 
    avg_tto = mean(tto_vals)
end
avg_ttb = mean(ttb_vals)

avg_profit = mean(profit_vals)

println("\n================ SUMMARY STATISTICS ================\n")

@printf("%-20s %-15s\n", "Metric", "Average")

println("-"^40)

@printf("%-20s %-15.4f\n", "Objective", avg_obj)
@printf("%-20s %-15.4f\n", "Profit", avg_profit)
@printf("%-20s %-15.2f\n", "Iterations", avg_iter)
@printf("%-20s %-15.4f\n", "Time to Gap", avg_ttg)
@printf("%-20s %-15.4f\n", "Time to Optimal", avg_tto)
@printf("%-20s %-15.4f\n", "Time to Best", avg_ttb)