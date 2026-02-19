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
using JuMP, Gurobi, LinearAlgebra, Graphs, Plots
using JSON, PyCall, GeoInterface, LibGEOS, CSV, DataFrames, Statistics

# ---- Get data ---- #

Data_file_name = "Denmark_FictiveData_LonLat.jl" # Change this to "BiggerData.jl" for the larger dataset

data_file = get(ENV, "MODEL_DATA", joinpath(@__DIR__, "Data", Data_file_name))
if isfile(data_file)
    include(data_file)   # defines airport_coords, Population_coords, Pi, tau, etc.
else
    error("Data file not found: $data_file")
end

Max_Path_Stops = 4 # Maximum number of stops allowed in a path (excluding start and end), used for filtering reasonable paths

# --- Geographic distance helpers --- #
const EARTH_RADIUS_KM = 6371.0
deg2rad(deg) = deg * (pi/180)
function haversine(lon1, lat1, lon2, lat2)
    φ1 = deg2rad(lat1); φ2 = deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1); Δλ = deg2rad(lon2 - lon1)
    a = sin(Δφ/2)^2 + cos(φ1)*cos(φ2)*sin(Δλ/2)^2
    return 2 * EARTH_RADIUS_KM * asin(min(1, sqrt(a)))
end

# ---- Airport Data ---- #
N_i = Dict(i => String[] for i in keys(airport_coords))
N = collect(keys(airport_coords))
num_airports = length(N)


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
Dist_to_airports = Dict(k => Dict(i => begin
            ax, ay = airport_coords[i]
            px, py = Population_coords[k]
            # assume stored as (lon, lat)
            haversine(ax, ay, px, py)
        end for i in N) for k in keys(Population_coords))
                               
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
        a_lon, a_lat = airport_coords[rev_idx_airport[i]]
        b_lon, b_lat = airport_coords[rev_idx_airport[j]]
        d = haversine(a_lon, a_lat, b_lon, b_lat)
        
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

    # ----- Added filtering to remove paths that are too long (more than 3 nodes), unless they are the only option ----- #
    # Remove paths with many nodes, unless they're the only option for a (k,d) pair
    if !isempty(all_paths)
        shortest_idx = argmin(length.(all_paths))
        shortest_path = all_paths[shortest_idx]
    end
    
    filter!(p -> length(p) <= Max_Path_Stops, all_paths)
    if isempty(all_paths) && !isempty(all_paths)  # If all paths were filtered out
        all_paths = [shortest_path]
    end
    # ----- End of added filtering ----- #

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

model = Model(Gurobi.Optimizer)

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
# @objective(model, Max, mu[1]*sum(Pi[k]*phi[idx_Population[k], idx_destination[d]] for k in K for d in D))
# @objective(model, Min, mu[2]*sum(c[i]*y[idx_airport[i]] for i in N))

@objective(model, Max, mu[1]*sum(Pi[k]*phi[idx_Population[k], idx_destination[d]] for k in K for d in D) - 
                    mu[2]*sum(c[i]*y[idx_airport[i]] for i in N))


# ---- Constraints ---- #
    # @constraint(model, [i in N], rho[idx_airport[i]] == SP(y[idx_airport[i]], rho, i))

@constraint(model, [e in E], dist_matrix[idx_airport[e[1]], idx_airport[e[2]]] + rho[idx_airport[e[1]]] + rho[idx_airport[e[2]]] - M2 * (1 - z[idx_airport[e[1]], idx_airport[e[2]]]) <= tau)

@constraint(model, [p in P, e in E_p[p]], psi[idx_paths[p]] <= z[idx_airport[e[1]], idx_airport[e[2]]])

for k in K, d in D
    if d in Closest_airports[k]
        continue
    end
    @constraint(model, phi[idx_Population[k], idx_destination[d]] <= sum(psi[idx_paths[p]] for p in P_k_d[(k, d)]))
end

# @constraint(model, [k in K, d in D], phi[idx_Population[k], idx_destination[d]] <= sum(psi[idx_paths[p]] for p in P_k_d[(k, d)]))

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
    if length(E) <= 10
        println("z (Feasible edges) = ", value.(z))
    else
        println("z (Feasible edges) = [Too many edges to display]")
    end
    println("rho (Shortest dist to charger) = ", value.(rho))
    if length(P) <= 10
        println("psi (Path feasibility) = ", value.(psi))
    else
        println("psi (Path feasibility) = [Too many paths to display]")
    end
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
    node_marker = [rev_idx_airport[i] in D ? :square : :circle for i in 1:num_airports]
    node_size  = 8

    pop_coords = [Population_coords[k] for k in K]
    pop_xs = first.(pop_coords); pop_ys = last.(pop_coords)
    pop_labels = [k for k in K]
    
    plt = scatter(xs, ys, markersize=node_size, markercolor=node_color, markerstrokewidth=0, legend=false)

    # Plot non-destination airports
    other_idx = [i for i in 1:num_airports if !(rev_idx_airport[i] in D)]
    if !isempty(other_idx)
        scatter!(plt, xs[other_idx], ys[other_idx], markersize=node_size, markercolor=node_color[other_idx], marker=:circle, markerstrokewidth=0, legend=false)
    end

    # Plot destination airports with a star marker, larger and distinct color
    dest_idx = [i for i in 1:num_airports if rev_idx_airport[i] in D]
    if !isempty(dest_idx)
        scatter!(plt, xs[dest_idx], ys[dest_idx], markersize=node_size, markercolor=node_color[dest_idx], marker=:square, markerstrokewidth=0, legend=false)
    end

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
        col = :gray
        plot!(plt, [xs[i], xs[j]], [ys[i], ys[j]], color=col, linewidth=1.5)
    end

    for p in P
        if value(psi[idx_paths[p]]) > 0.5
            for e in E_p[p]
                i, j = idx_airport[e[1]], idx_airport[e[2]]
                plot!(plt, [xs[i], xs[j]], [ys[i], ys[j]], color=:blue, linewidth=2)
            end
        end
    end

    for i in 1:num_airports
        annotate!(xs[i]+5, ys[i]+5, rev_idx_airport[i])
    end
    for k in K
        pop_x, pop_y = pop_coords[idx_Population[k]]
        annotate!(pop_x+5, pop_y+5, pop_labels[idx_Population[k]])
    end

    # Ensure output directory exists and save plot
    outdir = joinpath(@__DIR__, "Solution Plots")
    isdir(outdir) || mkpath(outdir)
    savefig(plt, joinpath(outdir, "solution_plot_$(Data_file_name).png"))
end

plot_simple()


## ---- Folium Map Visualization ---- #
include(joinpath(@__DIR__, "MapPlot.jl"))

routes = Vector{Vector{Tuple{Float64,Float64}}}()
for p in P
    if value(psi[idx_paths[p]]) > 0.5
        push!(routes, [airport_coords[rev_idx_airport[i]] for i in p])
    end
end

charging_airports = [rev_idx_airport[i] for i in 1:num_airports if value(y[i]) > 0.5]
population_links = Vector{Tuple{String,String}}()
for k in K
    for airport in Closest_airports[k]
        push!(population_links, (k, airport))
    end
end

save_folium_map(airport_coords, Population_coords, routes, charging_airports, population_links, D)

