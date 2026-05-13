using XLSX
using DataFrames
using Random
using Distributions
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

excel_file = joinpath("inputData/inputDataGiant.xlsx")
parameter_file = joinpath("inputData/Parameters.xlsx")
data = load_data(excel_file, parameter_file)

Vmax = maximum(data.V)
rt = zeros(Vmax, Vmax)
rt = zeros(Vmax, Vmax)
for i in data.V, j in data.V
    rt[i,j] = data.rt[(i,j)]
    rt[i,j] = data.rt[(i,j)]
end


### --------------------------- ###

# Test target solution feasibility
println("\n========== TESTING TARGET SOLUTION FEASIBILITY ==========\n")

# Create the target solution
evtol1 = planeSolution(2, [1, 3, 1], [47, 34])
evtol2 = planeSolution(0, [1], Int64[])
evtol3 = planeSolution(3, [2, 5, 4, 2], [0, 30, 30])
evtol4 = planeSolution(2, [2, 3, 2], [0, 38])

target_solution = allPlaneSolution([evtol1, evtol2, evtol3, evtol4])

println("Target Solution:")
println("  eVTOL1: $(evtol1.flightLegs) | $(Int.(evtol1.route)) | $(Int.(evtol1.turnaroundTime))")
println("  eVTOL2: $(evtol2.flightLegs) | $(Int.(evtol2.route)) |")
println("  eVTOL3: $(evtol3.flightLegs) | $(Int.(evtol3.route)) | $(Int.(evtol3.turnaroundTime))")
println("  eVTOL4: $(evtol4.flightLegs) | $(Int.(evtol4.route)) | $(Int.(evtol4.turnaroundTime))")

# Test feasibility using FeasibilityCheck (P function)
println("\nRunning feasibility check...")
P = FeasibilityCheck(
    Float32(data.bmax),
    Float32(data.bmid),
    Float32(data.bmin),
    data.dist,
    Float32(data.ec),
    Float32(data.battery_per_km),
    target_solution,
    Int.(rt),
    Int(round(data.ET)),
    Int(maximum(data.T)),
    Int(maximum(data.V)),
    data.cap_v,
    data.b_penalty
)

println("\nPenalty Vector P:")
println("  P = $P")
println("  Sum(P) = $(sum(P))")

if sum(P) < 1_000_000
    println("✓ Solution is FEASIBLE (penalty < 1,000,000)")
else
    println("✗ Solution is INFEASIBLE (penalty ≥ 1,000,000)")
end

println("\n========================================================\n")

Objective = obj(target_solution, data, Int.(rt))

println("Objective =  $(Objective)")