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

excel_file     = joinpath(@__DIR__, "..", "..", "inputData", "Experiments/inputDataEx10_15_35.xlsx")
parameter_file = joinpath(@__DIR__, "..", "..", "inputData", "Parameters.xlsx")

data = load_data(excel_file, parameter_file, excel_file, 60)

Vmax = maximum(data.V)
rt = zeros(Vmax, Vmax)
rt = zeros(Vmax, Vmax)
for i in data.V, j in data.V
    rt[i,j] = data.rt[(i,j)]
    rt[i,j] = data.rt[(i,j)]
end

maxTurnaround = Int64(data.ET)
Maxtime = Int32(100)
top_c = 4
Optimum = Inf
# Random.seed!(1234)

LoopIter = 5

obj_vals = zeros(LoopIter)
iter_vals = zeros(LoopIter)
ttg_vals = zeros(LoopIter)
tto_vals = zeros(LoopIter)
ttb_vals = zeros(LoopIter)
profit_vals = zeros(LoopIter)
Passenger_vals = zeros(LoopIter)

for i in 1:LoopIter
    (obj_vals[i], best_sols, iter_vals[i], ttg_vals[i], tto_vals[i], ttb_vals[i], profit_vals[i]) =
        HeuristicSA(maxTurnaround, Maxtime, data, rt, top_c; optimal_obj = Float32(Optimum))
    Passenger_vals[i] = passenger_demand_served(best_sols, data, rt)
    assignments, scheduled = assign_passengersV2(best_sols, data, Int.(rt))
    println("Best solution:")
    println("Objective Value: $(obj_vals[i])")
    println("Profit Value: $(profit_vals[i])")
    print_chromosome_table(best_sols)
    print_assignments(assignments, data)
end


avg_obj = mean(obj_vals)
avg_iter = mean(iter_vals)
avg_passenger = mean(Passenger_vals)
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

std_obj = std(obj_vals)
std_iter = std(iter_vals)
std_passenger = std(Passenger_vals)

if Inf in ttg_vals
    std_ttg = Inf32
else
    std_ttg = std(ttg_vals)
end

if Inf in tto_vals
    std_tto = Inf32
else
    std_tto = std(tto_vals)
end

std_ttb = std(ttb_vals)
std_profit = std(profit_vals)

println("\n================ SUMMARY STATISTICS ================\n")

@printf("%-20s %-15s %-15s\n", "Metric", "Average", "Std Dev")

println("-"^55)

@printf("%-20s %-15.4f %-15.4f\n", "Objective", avg_obj, std_obj)
@printf("%-20s %-15.4f %-15.4f\n", "Profit", avg_profit, std_profit)
@printf("%-20s %-15.2f %-15.2f\n", "Iterations", avg_iter, std_iter)
@printf("%-20s %-15.4f %-15.4f\n", "Time to Gap", avg_ttg, std_ttg)
@printf("%-20s %-15.4f %-15.4f\n", "Time to Optimal", avg_tto, std_tto)
@printf("%-20s %-15.4f %-15.4f\n", "Time to Best", avg_ttb, std_ttb)
@printf("%-20s %-15.4f %-15.4f\n", "Served Passengers", avg_passenger, std_passenger)
