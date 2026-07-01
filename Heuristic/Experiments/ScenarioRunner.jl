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

# ── Scenario table ────────────────────────────────────────────────────────────
# (scenario_name, parameter_file, drive_radius_minutes)
const SCENARIOS = [
    ("Base", "Parameters Base.xlsx",  60),
    ("s1",   "Parameters s1.xlsx",    60),
    ("s2",   "Parameters s2.xlsx",    60),
    ("s3",   "Parameters s3.xlsx",   120),
    ("s4",   "Parameters s4.xlsx",    60),
    ("s5",   "Parameters s5.xlsx",    60),
    ("s6",   "Parameters s6.xlsx",   120),
    ("s7",   "Parameters s7.xlsx",   120),
    ("s8",   "Parameters s8.xlsx",   120),
    ("s9",   "Parameters s9.xlsx",   120),
    ("s10",  "Parameters s10.xlsx",  120),
]

const SCENARIO_PARAMS_DIR = joinpath(@__DIR__, "..", "..", "inputData", "Scenario Parameters")
const INFRA_FILE          = joinpath(@__DIR__, "..", "..", "inputData", "inputDataHumongous.xlsx")
const DEMAND_FILE         = joinpath(@__DIR__, "..", "..", "inputData", "LTM_demand.xlsx")

# ── Helpers ───────────────────────────────────────────────────────────────────

function get_et(parameter_file::String)
    et = nothing
    XLSX.openxlsx(parameter_file) do xf
        for row in XLSX.eachrow(xf["Parameters"])
            k, v = row[1], row[2]
            if k !== missing && v !== missing && String(k) == "ET"
                et = Int64(round(parse(Float64, replace(string(v), "," => "."))))
                break
            end
        end
    end
    et === nothing && error("ET not found in $parameter_file")
    return et
end

function load_scenario_data(param_file::String, drive_radius::Int)
    run(`python c:/Users/Kapta/Documents/DTU/Speciale/GitHubFiles/DataGen/LTM/LTMPassAssignRand.py`)
    data = load_data(INFRA_FILE, param_file, DEMAND_FILE, drive_radius)

    vmax = maximum(data.V)
    rt = zeros(vmax, vmax)
    for i in data.V, j in data.V
        rt[i, j] = data.rt[(i, j)]
    end
    return data, rt
end

# ── Batch runner for one scenario ─────────────────────────────────────────────

function run_scenario_batch(scenario_name::String, param_file::String, drive_radius::Int;
                            n_runs::Int, maxtime::Int32, top_c::Int, seed::Int)
    max_turnaround = get_et(param_file)

    best_objs   = Float64[]
    best_sols   = allPlaneSolution[]
    data_list   = Any[]
    rt_list     = Any[]
    iter_list   = Int[]
    profit_list = Float64[]

    for run_i in 1:n_runs
        Random.seed!(seed + run_i - 1)
        data, rt = load_scenario_data(param_file, drive_radius)

        best_obj, best_sol, iters, _, _, _, profit =
            HeuristicSA(max_turnaround, maxtime, data, rt, top_c)

        push!(best_objs,   best_obj)
        push!(best_sols,   best_sol)
        push!(data_list,   data)
        push!(rt_list,     rt)
        push!(iter_list,   iters)
        push!(profit_list, profit)

        @printf("  [%s] Run %2d/%d  obj=%8.1f  profit=%8.1f  iters=%d\n",
                scenario_name, run_i, n_runs,
                best_obj, profit, iters)
    end

    return best_objs, best_sols, data_list, rt_list, iter_list, profit_list
end

# ── Save results for one scenario ─────────────────────────────────────────────

function save_scenario_results(scenario_name::String,
                               best_objs, best_sols, data_list, rt_list,
                               iter_list, profit_list;
                               out_dir::String)
    mkpath(out_dir)

    # Objective + profit per run
    CSV.write(
        joinpath(out_dir, "runs.csv"),
        DataFrame(run        = 1:length(best_objs),
                  best_obj   = best_objs,
                  profit     = profit_list,
                  iterations = iter_list)
    )

    # KPIs for every run
    kpi_rows = [
        merge(
            (scenario    = scenario_name,
             run         = i,
             best_obj    = best_objs[i],
             best_profit = profit_list[i]),
            solution_kpis(best_sols[i], data_list[i], rt_list[i])
        )
        for i in eachindex(best_sols)
    ]
    CSV.write(joinpath(out_dir, "kpis.csv"), DataFrame(kpi_rows))

    # Best solution (highest profit run) written to a text file
    best_idx = argmax(profit_list)
    open(joinpath(out_dir, "best_solution.txt"), "w") do io
        println(io, "Scenario  : $scenario_name")
        println(io, "Best run  : $best_idx / $(length(best_objs))")
        println(io, "Obj value : $(best_objs[best_idx])")
        println(io, "Profit    : $(profit_list[best_idx])")
        println(io)

        redirect_stdout(io) do
            print_chromosome_table(best_sols[best_idx])
            println()
            assignments, scheduled =
                assign_passengersV2(best_sols[best_idx],
                                    data_list[best_idx],
                                    Int.(rt_list[best_idx]))
            print_assignments(assignments, data_list[best_idx])
            println()
            print_schedule_pretty(scheduled)
        end
    end

    println("  Saved: $out_dir")
    return kpi_rows
end

# ── Top-level entry point ─────────────────────────────────────────────────────

function run_all_scenarios(;
    n_runs  :: Int   = 20,
    maxtime :: Int32 = Int32(200),
    top_c   :: Int   = 4,
    seed    :: Int   = 1,
)
    results_root = joinpath(@__DIR__, "Results", "Scenarios")
    mkpath(results_root)

    all_kpi_rows = Any[]

    for (scenario_name, param_filename, drive_radius) in SCENARIOS
        println("\n═══════════════════════════════════════════════════════")
        println("  Scenario: $scenario_name   Drive_radius: $drive_radius min")
        println("═══════════════════════════════════════════════════════")

        param_file = joinpath(SCENARIO_PARAMS_DIR, param_filename)

        best_objs, best_sols, data_list, rt_list, iter_list, profit_list =
            run_scenario_batch(scenario_name, param_file, drive_radius;
                               n_runs=n_runs, maxtime=maxtime,
                               top_c=top_c, seed=seed)

        out_dir = joinpath(results_root, scenario_name)
        kpi_rows = save_scenario_results(scenario_name,
                                         best_objs, best_sols,
                                         data_list, rt_list,
                                         iter_list, profit_list;
                                         out_dir=out_dir)
        append!(all_kpi_rows, kpi_rows)

        best_profit = maximum(profit_list)
        @printf("  Done.  max_profit=%.1f  median_obj=%.1f\n",
                best_profit, median(best_objs))
    end

    # Combined KPI table across all scenarios
    combined_csv = joinpath(results_root, "all_scenarios_kpis.csv")
    CSV.write(combined_csv, DataFrame(all_kpi_rows))
    println("\nAll scenarios complete.")
    println("Combined KPIs: $combined_csv")
end

run_all_scenarios()
