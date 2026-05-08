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

candiateroutes = Candidate_Route(data)

iter = 10
objs = zeros(iter)
P_values = []
evtols_init_values = []
results = []




for i in 1:iter

    evtols_init = Construction_Heuristic(data, candiateroutes)
    push!(evtols_init_values, evtols_init)

    assignments, scheduled = assign_passengersV2(evtols_init, data, Int.(rt))

    P = FeasibilityCheck(Float32(data.bmax),Float32(data.bmid), Float32(data.bmin),data.dist,Float32(data.ec),Float32(data.battery_per_km),
                evtols_init,Int.(rt),Int(round(data.ET)),maximum(Int.(data.T)),maximum(data.V),data.cap_v, data.b_penalty)

    fitness = fitnessFunction(evtols_init,assignments,Float32(data.bmax), Float32(data.bmid), Float32(data.bmin),data.dist, Float32(data.ec),
                    Float32(data.battery_per_km), Int.(rt), Int(round(data.ET)), maximum(Int.(data.T)), maximum(data.V), data.cap_v, data)

    objs[i] = fitness
    push!(P_values, P)
    push!(results, (fitness=fitness, P=P, evtols_init=evtols_init, assignments=assignments, scheduled=scheduled))
end

sorted_indices = sortperm(objs, rev=true)
println("Top 10 entries:")
println(objs[sorted_indices[1:10]])

sorted_results = sort(results, by = x -> x.fitness, rev = true)

println("\nTop 5 entries with associated P values:")
for i in 1:5
    r = sorted_results[i]
    println("Fitness: $(r.fitness), P: $(r.P)")
end

zero_p_count = sum(sum(P) == 0 for P in P_values)
zero_p_percentage = (zero_p_count / length(P_values)) * 100
println("\nPercentage of objs with sum(P) == 0: $(zero_p_percentage)%")

# For solutions with P != 0, print which indices of P are nonzero and overall distribution
nonzero_indices_list = []
max_len = maximum(length.(P_values))
counts = zeros(Int, max_len)
for (run_idx, P) in enumerate(P_values)
    if sum(P) != 0
        inds = findall(!=(0), P)
        push!(nonzero_indices_list, (run=run_idx, indices=inds))
        for i in inds
            counts[i] += 1
        end
        # println("Run $(run_idx): nonzero P indices -> $(inds)")
    end
end

println("\nDistribution (count of runs where P[index] != 0):")
for i in 1:length(counts)
    if counts[i] > 0
        println("Index $(i): $(counts[i])")
    end
end


println("hey")