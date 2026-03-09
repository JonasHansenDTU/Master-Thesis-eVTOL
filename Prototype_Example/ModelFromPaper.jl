"""
Table 1 — Model notation.

Inputs: sets & parameters

- N : set of airport nodes, indexed by i and j
- N_i : set of airport nodes adjacent to node i, i.e., j ∈ N such that d_ij ≤ τ and i ≠ j
- E : set of flight edges, indexed by (i, j)
- K : set of population areas, indexed by k
- P : set of paths, indexed by p
- D : set of destinations, indexed by d
- P_k^d : subset of paths from population area k to destination d
- E_p : subset of flight edges in path p
- d_ij : travel distance between i and j
- τ : maximum travel range on a single charge
- c_i : cost of activating node i as a charging base
- π_k : population living in k
- μ_1, μ_2 : weights of the two objectives in the objective function

Variables

- y_i ∈ {0,1} : = 1 if node i is activated as a charging base
- z_ij ∈ {0,1} : = 1 if flight edge (i, j) is feasible
- ρ_i ∈ R^+ : value of the shortest path from i to a charging base
- ψ_p ∈ {0,1} : = 1 if path p can be used
- φ^d_k ∈ {0,1} : = 1 if area k is covered w.r.t. destination d

"""
# ---- Settings/Fixed values ---- #
using JuMP, HiGHS, LinearAlgebra, Graphs, Plots

tau = 400 
Reasnable_range_factor = 2
Travel_range_to_airort_factor = 1.74 # How much farther than the closest airport the population is willing to travel
mu = [1.0, 0.5] # Weights for the objective function (coverage vs cost)

M1 = 1000 # A sufficiently large constant for Big-M constraints
M2 = 1000 # Another large constant for linearization of the SP function

# ---- Airport Data ---- #
airport_coords = Dict("A"=>(230,140), "B"=>(300,250), "C"=>(430,250), "D"=>(620,160))
N_i = Dict(i => String[] for i in keys(airport_coords))
N = collect(keys(airport_coords))
num_airports = length(N)
D = ["D"] # Destination set
c = Dict("A"=>10, "B"=>20, "C"=>15, "D"=>25) # Cost of activating each airport as a charging base

idx_airport = Dict(name => i for (i, name) in enumerate(N))
rev_idx_airport = Dict(i => name for (i, name) in enumerate(N))

idx_destination = Dict(name => i for (i, name) in enumerate(D))
rev_idx_destination = Dict(i => name for (i, name) in enumerate(D))

# ---- Edges and Paths ---- #
E = []
E_p = Dict{Vector{Int}, Vector{Tuple{String,String}}}() # For each path (vector of nodes), store the edges (tuple of node names) in that path
P = []
g = SimpleDiGraph(num_airports)

# ---- Population Data ---- #
Population_coords = Dict("K1"=>(60,60), "K2"=>(80,330))
Pi = Dict("K1"=>1000, "K2"=>1500) # Population in each area
Dist_to_airports = Dict(k => Dict(i => hypot(airport_coords[i][1]-Population_coords[k][1], 
                                  airport_coords[i][2]-Population_coords[k][2]) for i in N) for k in keys(Population_coords))
                               
idx_Population = Dict(name => i for (i, name) in enumerate(keys(Population_coords)))
rev_idx_Population = Dict(i => name for (i, name) in enumerate(keys(Population_coords)))

K = collect(keys(Population_coords))
Closest_airports = Dict(k => begin
    distances = Dist_to_airports[k]
    min_dist = minimum(values(distances))
    filter_dist = min_dist * Travel_range_to_airort_factor
    [airport for airport in N if distances[airport] <= filter_dist]
end for k in K) # For each population area, find the closest airports

P_k_d = Dict{Tuple{String,String}, Vector{Vector{Int}}}() # For each (k, d), store the reasonable paths from k to d

# ---- Computaion of Distrance Matrix


# Find all feasible EDGES based on Range (tau) 
dist_matrix = zeros(num_airports, num_airports)
for i in 1:num_airports, j in 1:num_airports
    if i != j
        d = hypot(airport_coords[rev_idx_airport[i]][1]-airport_coords[rev_idx_airport[j]][1], 
                  airport_coords[rev_idx_airport[i]][2]-airport_coords[rev_idx_airport[j]][2])
        
        dist_matrix[i, j] = d
        if d <= tau
            push!(N_i[rev_idx_airport[i]], rev_idx_airport[j])
            push!(E, (rev_idx_airport[i], rev_idx_airport[j]))
            add_edge!(g, i, j)
        end
    end
end

# ---- Generate reasonable paths for each (population area, destination) pair ---- #

function get_all_reasonable_paths(start_node, end_node, max_d)
    # Using a simple DFS to find paths under the travel limit
    all_paths = []
    function find_paths(curr, path, d)
        if curr == end_node
            push!(all_paths, copy(path))
            return
        end
        for nb in neighbors(g, curr)
            new_d = d + dist_matrix[curr, nb]
            if !(nb in path) && new_d <= max_d
                push!(path, nb)
                find_paths(nb, path, new_d)
                pop!(path)
            end
        end
    end
    find_paths(start_node, [start_node], 0.0)
    return all_paths
end


for d in D
    for i in N
        if i == d
            continue
        end
        reasnable_range = Reasnable_range_factor * dist_matrix[idx_airport[i], idx_airport[d]]
        k_paths = get_all_reasonable_paths(idx_airport[i], idx_airport[d], reasnable_range)
        append!(P, k_paths)
        for p in k_paths
            E_p[p] = [(rev_idx_airport[p[j]], rev_idx_airport[p[j+1]]) for j in 1:(length(p)-1)]
        end

        for k in K
            if i in Closest_airports[k]
                # append new paths to existing entry (or create it)
                append!(get!(P_k_d, (k, d), Vector{Vector{Int}}()), k_paths)
            end
        end
    end
end

idx_paths = Dict(p => i for (i, p) in enumerate(P))
rev_idx_paths = Dict(i => p for (i, p) in enumerate(P))

function SP(y, rho, i)
    sp = minimum(M1*(1-y), minimum(dist_matrix[idx_airport[i],idx_airport[j]]+rho[idx_airport[j]] for j in N_i[i]))
    return sp
end

## ---- 4. MODEL (EACN-REG Formulation) ---- #

model = Model(HiGHS.Optimizer)

# ---- Variables ---- #
@variable(model, y[1:num_airports], Bin) # Charging Bases
@variable(model, z[1:length(E), 1:length(E)], Bin) # Feasible edges
@variable(model, rho[1:num_airports] >= 0) # Shortest dist to charger
@variable(model, psi[1:length(P)], Bin) # Path feasibility
@variable(model, phi[1:length(K), 1:length(D)], Bin) # Population coverage

#Varaibles to linarize the SP function
@variable(model, w[1:num_airports, 1:num_airports], Bin) # w[i,j] = 1 if the shortest path for node i goes through node j
@variable(model, x[1:num_airports], Bin) # x[i] = 1 if no feasible path from node i to acharging station exists


# ---- Objectives ---- #
@objective(model, Max, mu[1]*sum(Pi[k]*phi[idx_Population[k], idx_destination[d]] for k in K for d in D))
# @objective(model, Min, mu[2]*sum(c[i]*y[idx_airport[i]] for i in N))



# ---- Constraints ---- #
    # @constraint(model, [i in N], rho[idx_airport[i]] == SP(y[idx_airport[i]], rho, i))

@constraint(model, [e in E], dist_matrix[idx_airport[e[1]], idx_airport[e[2]]] + rho[idx_airport[e[1]]] + rho[idx_airport[e[2]]] - M2 * (1 - z[idx_airport[e[1]], idx_airport[e[2]]]) <= tau)

@constraint(model, [p in P, e in E_p[p]], psi[idx_paths[p]] <= z[idx_airport[e[1]], idx_airport[e[2]]])

@constraint(model, [k in K, d in D], phi[idx_Population[k], idx_destination[d]] <= sum(psi[idx_paths[p]] for p in P_k_d[(k, d)]))

# ---- Shortest Path Constraints (Linearized) ---- #
@constraint(model, [i in N], rho[idx_airport[i]] <= M1 * (1 - y[idx_airport[i]])) # 
@constraint(model, [i in N, j in N_i[i]], rho[idx_airport[i]] 
                                            <= dist_matrix[idx_airport[i], idx_airport[j]] + 
                                            rho[idx_airport[j]])

@constraint(model, [i in N, j in N_i[i]], rho[idx_airport[i]] 
                                            >= dist_matrix[idx_airport[i], idx_airport[j]] + 
                                            rho[idx_airport[j]] - M2 * (1 - w[idx_airport[i], idx_airport[j]]))

@constraint(model, [i in N], rho[idx_airport[i]] >= M1 * x[idx_airport[i]]) # If x[i] = 1, then rho[i] must be at least M1 (i.e., no feasible path)
@constraint(model, [i in N], sum(w[idx_airport[i], idx_airport[j]] for j in N_i[i]) + 
                            y[idx_airport[i]] + x[idx_airport[i]] == 1) # Either one path is chosen or no path is feasible


## ---- Optimize and Display Results ---- #
optimize!(model)

if termination_status(model) == OPTIMAL
    println("Optimal solution found!")
    println("Objective value: ", objective_value(model))
    println("y (Charging Bases) = ", value.(y))
    println("z (Feasible edges) = ", value.(z))
    println("rho (Shortest dist to charger) = ", value.(rho))
    println("psi (Path feasibility) = ", value.(psi))
    println("phi (Population coverage) = ", value.(phi))
else
    println("No optimal solution found. Status: ", termination_status(model))
end


## ---- Visualization of results ---- #

using Plots

function plot_simple()
    coords = [airport_coords[rev_idx_airport[i]] for i in 1:num_airports]
    xs = first.(coords); ys = last.(coords)
    node_color = [value(y[i]) > 0.5 ? :red : :lightgray for i in 1:num_airports]
    node_size  = 8

    pop_coords = [Population_coords[k] for k in K]
    pop_xs = first.(pop_coords); pop_ys = last.(pop_coords)
    
    plt = scatter(xs, ys, markersize=node_size, markercolor=node_color, legend=false)
    scatter!(plt, pop_xs, pop_ys, markersize=6, markercolor=:green, markerstrokewidth=0, label="Population")

    for k in K
        pop_x, pop_y = pop_coords[idx_Population[k]]
        for airport in Closest_airports[k]
            airport_idx = idx_airport[airport]
            plot!(plt, [pop_x, xs[airport_idx]], [pop_y, ys[airport_idx]], color=:green, linestyle=:dash, linewidth=1)
        end
    end

    for e in collect(edges(g))
        i,j = src(e), dst(e)
        col = (value(z[i,j]) > 0.5) ? :blue : :gray
        plot!(plt, [xs[i], xs[j]], [ys[i], ys[j]], color=col, linewidth=1.5)
    end
    for i in 1:num_airports
        annotate!(xs[i]+5, ys[i]+5, rev_idx_airport[i])
    end
    savefig(plt, "solution_plot.png")
end

plot_simple()