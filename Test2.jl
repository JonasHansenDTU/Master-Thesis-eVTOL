using JuMP, HiGHS, LinearAlgebra, Graphs

# --- 1. DATA & SETTINGS ---
airport_coords = Dict("A"=>(0,0), "B"=>(100,0), "C"=>(200,0), "D"=>(100,50))
range_tau = 250.0
max_travel_dist = 400.0 # Paper's "Reasonable Travel" Filter [cite: 96]

BigM = 1000.0 # A large constant for Big-M constraints




# --- 2. AUTOMATIC PRE-PROCESSING ---
airports = collect(keys(airport_coords))
n = length(airports)
idx = Dict(name => i for (i, name) in enumerate(airports))
rev_idx = Dict(i => name for (i, name) in enumerate(airports))

# Step A: Find all feasible EDGES based on Range 
dist_matrix = zeros(n, n)
g = SimpleDiGraph(n)
for i in 1:n, j in 1:n
    if i != j
        d = hypot(airport_coords[rev_idx[i]][1]-airport_coords[rev_idx[j]][1], 
                  airport_coords[rev_idx[i]][2]-airport_coords[rev_idx[j]][2])
        if d <= range_tau
            dist_matrix[i, j] = d
            add_edge!(g, i, j)
        end
    end
end

# Step B: Find all reasonable PATHS (Automatic Routing) [cite: 435, 439]
# We find all paths from 'A' to 'C' that aren't too long.
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

# Generate paths for k1 (Location A) to reach Destination C
k1_paths = get_all_reasonable_paths(idx["A"], idx["C"], max_travel_dist)

# --- 3. MODEL (EACN-REG Formulation) ---
model = Model(HiGHS.Optimizer)
set_silent(model)

@variable(model, y[1:n], Bin)          # Charging Bases [cite: 98]
@variable(model, rho[1:n] >= 0)        # Shortest dist to charger [cite: 99]
@variable(model, z[1:n, 1:n], Bin)     # Edge feasibility [cite: 102]
@variable(model, psi[1:length(k1_paths)], Bin) # Path feasibility [cite: 103]
@variable(model, phi, Bin)             # Population Coverage [cite: 104]

# Shortest Path & Edge Feasibility Logic [cite: 116, 125]
for i in 1:n
    @constraint(model, rho[i] <= range_tau * (1 - y[i])) # If no charger, rho must be less than range_tau; if charger exists, rho must be 0
    for j in 1:n
        if dist_matrix[i,j] > 0
            # the closest charger to node i is at most the closed distance to j plus the closest charger to j
            @constraint(model, rho[i] <= dist_matrix[i,j] + rho[j])
            # rho_i + d_ij + rho_j <= tau [cite: 116, 138] ensures that if edge (i,j) is used, then the distance plus the "detours" to chargers must be within range_tau
            @constraint(model, dist_matrix[i,j] + rho[i] + rho[j] <= range_tau + BigM*(1 - z[i,j]))
        end
    end
end

# Path Logic: Path is usable only if all its edges are feasible [cite: 116, 138]
for (p_idx, path_nodes) in enumerate(k1_paths)
    for i in 1:(length(path_nodes)-1)
        u, v = path_nodes[i], path_nodes[i+1]
        @constraint(model, psi[p_idx] <= z[u, v])
    end
end

@constraint(model, phi <= sum(psi)) # Covered if any path works [cite: 116]

for i in 1:n
    # 1. Collect neighbor distances and compute a safe minimum
    neighbor_dists = [dist_matrix[i,j] for j in 1:n if dist_matrix[i,j] > 0]
    if isempty(neighbor_dists)
        min_dist_to_neighbor = range_tau
    else
        min_dist_to_neighbor = minimum(neighbor_dists)
    end
    # Tight Big-M per node (non-negative)
    M1_i = max(range_tau - min_dist_to_neighbor, 0.0)

    # 2. Add the Upper Bound [cite: 149]
    @constraint(model, rho[i] <= M1_i * (1 - y[i]))

    # 3. Add the Lower Bound (prevents rho=0 when y=0) [cite: 159]
    @constraint(model, rho[i] >= min(min_dist_to_neighbor, M1_i) * (1 - y[i]))
end

@objective(model, Max, 500*phi - 10*sum(y))
optimize!(model)

# --- 4. RESULTS ---
println("Automated Search found $(length(k1_paths)) possible routes.")
println("Optimal Chargers: ", [rev_idx[i] for i in 1:n if value(y[i]) > 0.5])

# --- Pretty-print variable values ---
using Printf

println("\n================ Variable values ================")

println("\n-- Chargers `y` (1=open, 0=closed) --")
for i in 1:n
    @printf(" %-3s : %d\n", rev_idx[i], Int(round(value(y[i]))))
end

println("\n-- Distance to nearest charger `rho` (units) --")
for i in 1:n
    @printf(" %-3s : %6.2f\n", rev_idx[i], value(rho[i]))
end

println("\n-- Edge feasibility `z` (1=feasible) --")
for i in 1:n, j in 1:n
    if i != j && value(z[i,j]) > 0.5
        @printf(" %s -> %s : 1\n", rev_idx[i], rev_idx[j])
    end
end

println("\n-- Path feasibility `psi` --")
for (p_idx, path_nodes) in enumerate(k1_paths)
    path_names = join([rev_idx[n] for n in path_nodes], " -> ")
    @printf(" P%02d [%s] : %d\n", p_idx, path_names, Int(round(value(psi[p_idx]))))
end

println("\n-- Coverage `phi` --")
@printf(" phi : %d\n", Int(round(value(phi))))

println("\n-- Distance matrix (only feasible edges shown) --")
for i in 1:n, j in 1:n
    if dist_matrix[i,j] > 0
        @printf(" %s -> %s : %6.1f\n", rev_idx[i], rev_idx[j], dist_matrix[i,j])
    end
end

println("\n-- Candidate paths (with total path distance) --")
for (p_idx, path_nodes) in enumerate(k1_paths)
    path_names = [rev_idx[n] for n in path_nodes]
    pd = sum(dist_matrix[path_nodes[k], path_nodes[k+1]] for k in 1:(length(path_nodes)-1))
    @printf(" P%02d: %s  (dist=%.1f)\n", p_idx, join(path_names, " -> "), pd)
end

println("\n===============================================\n")