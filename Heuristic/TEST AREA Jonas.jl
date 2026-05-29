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

excel_file = joinpath("inputData/LTM_demand.xlsx")
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



plane = Vector{planeSolution}(undef, 1)

plane[1] = planeSolution(
    Int32(8),
    Int32[2, 4, 2, 5, 2, 3, 1, 4, 2],
    Int32[0, 30, 30, 50, 30, 30, 30, 30]
)

single_sol = allPlaneSolution(plane)

UpdateTurnAroundTimes(single_sol, 1, 400, data)


fitness, _ = score_single_plane_solution(plane[1], data, rt)
refined_fitness, refined_sol = local_search_single_plane_solution(single_sol, fitness, 400, data, rt)

fitness, _ = score_single_plane_solution(refined_sol.planes[1], data, rt)