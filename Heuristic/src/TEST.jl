# =========================
# Template: test OR solution feasibility
# =========================

# Assumes these are already available:
# - data = load_data(...)
# - FeasibilityCheck(...)
# - planeSolution / allPlaneSolution structs

src_dir = joinpath(@__DIR__)
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

excel_file = joinpath("inputData/inputDataMini.xlsx")
parameter_file = joinpath("inputData/Parameters.xlsx")
data = load_data(excel_file, parameter_file)

function build_rt_matrix(data)
    Vmax = maximum(data.V)
    rt = zeros(Int, Vmax, Vmax)
    for i in data.V, j in data.V
        rt[i, j] = data.rt[(i, j)]
    end
    return rt
end

"""
routes_by_plane: Dict(plane_id => [v1, v2, ..., vk])   # nodes visited, including start/end
waits_by_plane : Dict(plane_id => [t1, t2, ..., t{k-1}]) # turnaround per leg
"""
function test_or_solution_feasibility(data, routes_by_plane::Dict{Int,Vector{Int}},
                                      waits_by_plane::Dict{Int,Vector{Int}})
    planes = planeSolution[]

    # Build planeSolution in the same order as data.N
    for n in data.N
        route_int = get(routes_by_plane, n, Int[])
        waits_int = get(waits_by_plane, n, Int[])

        if isempty(route_int)
            # Empty flight plan: keep plane parked at base
            base = data.bv[n]
            push!(planes, planeSolution(Int32(0), Int32[base], Int32[]))
            continue
        end

        # Basic consistency checks
        if length(route_int) < 1
            error("Plane $(n): route must have at least one node.")
        end
        if length(waits_int) != length(route_int) - 1
            error("Plane $(n): waits length must equal route length - 1.")
        end

        push!(planes, planeSolution(
            Int32(length(route_int) - 1),
            Int32.(route_int),
            Int32.(waits_int)
        ))
    end

    evtols = allPlaneSolution(planes)
    rt = build_rt_matrix(data)

    P = FeasibilityCheck(
        Float32(data.bmax),
        Float32(data.bmin),
        data.dist,
        Float32(data.ec),
        Float32(data.battery_per_km),
        evtols,
        rt,
        Int(round(data.ET)),
        maximum(Int.(data.T)),
        maximum(data.V),
        data.cap_v
    )

    println("Feasibility vector P = ", P)
    println("Battery feasible:          ", P[1] == 0)
    println("Completion time feasible:  ", P[2] == 0)
    println("Vertiport capacity feasible:", P[3] == 0)
    println("Overall feasible:          ", all(P .== 0))

    return P, evtols
end

# =========================
# Example usage (replace with your OR solution)
# =========================
routes_by_plane = Dict(
    1 => [1, 3, 1],
    2 => [2],            # no flights (parked)
    3 => [5, 2, 3, 5]
)

waits_by_plane = Dict(
    1 => [2, 2],
    2 => Int[],
    3 => [3, 2, 4]
)

P, evtols = test_or_solution_feasibility(data, routes_by_plane, waits_by_plane)