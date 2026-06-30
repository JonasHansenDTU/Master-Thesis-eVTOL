using XLSX
using DataFrames
using Random
using JuMP
using Gurobi
using CSV
using MathOptInterface
using Printf
using Statistics
using PyCall

const MOI = MathOptInterface

src_dir = joinpath(@__DIR__, "..", "src")
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

function load_heuristic_data()
    run(`python c:/Users/Kapta/Documents/DTU/Speciale/GitHubFiles/DataGen/LTM/LTMPassAssignRand.py`)
    excel_file = joinpath(@__DIR__, "..", "..", "inputData", "inputDataHumongous.xlsx")
    parameter_file = joinpath(@__DIR__, "..", "..",  "inputData", "Parameters.xlsx")
    data = load_data(excel_file, parameter_file)

    vmax = maximum(data.V)
    rt = zeros(vmax, vmax)
    for i in data.V, j in data.V
        rt[i, j] = data.rt[(i, j)]
    end

    return data, rt
end

function load_turnaround_period()
    parameter_file = joinpath(@__DIR__, "..", "..", "inputData", "Parameters.xlsx")
    et = nothing

    XLSX.openxlsx(parameter_file) do xf
        sheet = xf["Parameters"]
        for row in XLSX.eachrow(sheet)
            key = row[1]
            val = row[2]
            if key !== missing && val !== missing && String(key) == "ET"
                et = Int64(round(parse(Float64, replace(string(val), "," => ".")))) 
                break
            end
        end
    end

    et === nothing && error("Could not find ET in Parameters.xlsx")
    return et
end

function run_heuristic_batch(n_runs::Int; max_turnaround::Int64, maxtime::Int32, top_c::Int, seed::Union{Nothing, Int}=nothing)
    best_objs = Float64[]
    best_sols = allPlaneSolution[]
    run_data_list = Any[]
    run_rt_list = Any[]
    iterations = Int[]
    profits = Float64[]

    for run in 1:n_runs
        if seed !== nothing
            Random.seed!(seed + run - 1)
        end

        data, rt = load_heuristic_data()
        best_obj, best_sol, iter_count, _, _, _, profit = HeuristicSA(max_turnaround, maxtime, data, rt, top_c)
        push!(best_objs, best_obj)
        push!(best_sols, best_sol)
        push!(run_data_list, data)
        push!(run_rt_list, rt)
        push!(iterations, iter_count)
        push!(profits, profit)


        println("Run $run / $n_runs -> best_obj = $(round(best_obj; digits=2)), iterations = $iter_count")
    end

    return best_objs, best_sols, run_data_list, run_rt_list, iterations, profits
end

function save_kpi_results(kpi_rows; out_csv::AbstractString)
    df = DataFrame(kpi_rows)
    CSV.write(out_csv, df)
    return df
end

function save_batch_results(best_objs::Vector{Float64}; out_csv::AbstractString)
    df = DataFrame(run = 1:length(best_objs), best_obj = best_objs)
    CSV.write(out_csv, df)
    return df
end

function save_boxplot(best_objs::Vector{Float64}; out_png::AbstractString, title_text::AbstractString)
    try
        pyimport("matplotlib.pyplot")
    catch
        pyimport_conda("matplotlib.pyplot", "matplotlib")
    end

    plt = pyimport("matplotlib.pyplot")

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.boxplot(best_objs, vert=true, patch_artist=true, labels=["best_obj"])
    ax.set_title(title_text)
    ax.set_ylabel("Objective value")
    ax.grid(true, axis="y", linestyle="--", alpha=0.35)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
end
# 
function main()
    max_turnaround = load_turnaround_period()
    maxtime = Int32(100)
    top_c = 4
    n_runs = 1

    out_dir = joinpath(@__DIR__, "Results")
    mkpath(out_dir)

    best_objs, best_sols, run_data_list, run_rt_list, iterations, profitvals = run_heuristic_batch(
        n_runs;
        max_turnaround=max_turnaround,
        maxtime=maxtime,
        top_c=top_c,
        seed=1,
    )

    best_idx = argmax(profitvals)
    best_sol = best_sols[best_idx]

    results_csv = joinpath(out_dir, "best_obj_runs.csv")
    save_batch_results(best_objs; out_csv=results_csv)

    results_csv = joinpath(out_dir, "best_profit_runs.csv")
    save_batch_results(profitvals; out_csv=results_csv)

    boxplot_png = joinpath(out_dir, "best_obj_boxplot.png")
    save_boxplot(
        best_objs;
        out_png=boxplot_png,
        title_text="Heuristic best_obj over $(n_runs) runs",
    )

    kpi_rows = [merge((run = i, best_obj = best_objs[i], best_profit = profitvals[i]), solution_kpis(best_sols[i], run_data_list[i], run_rt_list[i])) for i in eachindex(best_sols)]
    kpi_csv = joinpath(out_dir, "best_solution_kpis.csv")
    save_kpi_results(kpi_rows; out_csv=kpi_csv)

    println()
    println("Saved run results to: $results_csv")
    println("Saved boxplot to: $boxplot_png")
    println("Saved KPI results to: $kpi_csv")
    println("Summary: min=$(minimum(best_objs)), median=$(median(best_objs)), mean=$(mean(best_objs)), max=$(maximum(best_objs))")
    println("Iterations per run: $(iterations)\n")


    println("Best solution (higest profit):")
    println("Objective Value: $(best_objs[best_idx])")
    println("Profit Value: $(profitvals[best_idx])")
    print_chromosome_table(best_sol)

    assignments, scheduled = assign_passengersV2(best_sol, run_data_list[best_idx], Int.(run_rt_list[best_idx]))

    # println("Passenger Assignment")
    print_assignments(assignments, run_data_list[best_idx])

    # Export solution snapshots for visualization
    # snapshots = export_solution_snapshots(best_sol, scheduled, assignments, battery_levels, data)
    print_schedule_pretty(scheduled) 

end


main()
