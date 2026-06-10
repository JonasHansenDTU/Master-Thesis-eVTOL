###############################################################################
# RunMeasures.jl
#
# Standalone driver for the stochastic-programming measures VSS and EVPI.
#
#   RP   = E[profit] of the stochastic (two-stage) solution
#   EEV  = expected result of the expected-value (mean-scenario) solution
#   WS   = wait-and-see (per-scenario perfect-foresight) value
#   VSS  = RP - EEV     EVPI = WS - RP
#
# This is kept SEPARATE from MainCall_stochastic.jl so the plain two-stage run
# stays fast: the wait-and-see step solves every scenario deterministically and
# is the expensive part. Run this file only when you actually want the measures.
#
# It runs the two-stage heuristic ONCE to obtain RP, then computes EEV and WS.
# If you already have an RP you trust (e.g. the mean over a batch of runs), set
# USE_FIXED_RP = true and RP_FIXED below to skip the single two-stage solve.
###############################################################################

using XLSX
using DataFrames
using Random
using CSV
using Printf

src_dir = joinpath(@__DIR__, "..", "Heuristic", "src")

source_files = [
    "DataLoadFunc.jl",
    "FeasibilityFunc.jl",
    "InitialSolFunc.jl",
    "PassAssignFunc.jl",
    "FitnessFunc.jl",
    "NeighborhoodFunc.jl",
    "SANeighborhood.jl",
    "HeuristicFunc.jl",
]
for file in source_files
    include(joinpath(src_dir, file))
end

include(joinpath(@__DIR__, "scenario_generation_season.jl"))
include(joinpath(@__DIR__, "StochasticHeuristic.jl"))
include(joinpath(@__DIR__, "compute_vss_evpi.jl"))

###############################################################################
# Configuration (keep identical to MainCall_stochastic.jl for comparability)
###############################################################################

maxTurnaround = 100

MaxTime_1st        = Int32(90)
n_restarts         = 2
n_outer_iters      = 4
MaxTime_2nd_search = Int32(8)
MaxTime_2nd_final  = Int32(30)
price_boost        = 8.0
hard_penalty       = 50_000.0
top_c              = 4

# Budget for the deterministic per-scenario solves used by WS and EEV. Generous
# is better here so WS/EEV are well approximated and the ordering EEV<=RP<=WS
# is not violated by under-solving. Lower this (e.g. Int32(45)) to save time.
MaxTime_det        = Int32(90)

# If you already have a trustworthy RP (e.g. the mean over a batch of two-stage
# runs), set this to true and fill in RP_FIXED to skip the single solve below.
USE_FIXED_RP = false
RP_FIXED     = 0.0

###############################################################################
# Load data and generate scenarios (mirrors MainCall_stochastic.jl)
###############################################################################

println("Loading data …")
excel_file     = joinpath(@__DIR__, "..", "inputData", "inputDataGiant.xlsx")
parameter_file = joinpath(@__DIR__, "..", "inputData", "Parameters.xlsx")

data = load_data(excel_file, parameter_file)
time_per_km = derive_time_per_km(data)

rt_s, e_s, S, pi_s = generate_scenarios(
    data.V, data.lat, data.lon, data.rt, data.e, time_per_km
)

# Drop zero-probability scenarios for the active season (same as MainCall).
S_all    = collect(S)
S_active = [sc for sc in S_all if pi_s[sc] > 0.0]
pi_s     = Dict(sc => pi_s[sc] for sc in S_active)
S        = S_active
println("Scenarios (active season): $(length(S)) of $(length(S_all))")

###############################################################################
# Obtain RP
###############################################################################

if USE_FIXED_RP
    RP = RP_FIXED
    println("\nUsing fixed RP = $(round(RP, digits=2))")
else
    println("\nRunning the two-stage heuristic once to obtain RP …")
    result = stochastic_heuristic(
        data, rt_s, e_s, S, pi_s;
        maxTurnaround        = maxTurnaround,
        MaxTime_1st          = MaxTime_1st,
        MaxTime_2nd_search   = MaxTime_2nd_search,
        MaxTime_2nd_final    = MaxTime_2nd_final,
        top_c                = top_c,
        price_boost          = price_boost,
        hard_penalty         = hard_penalty,
        n_restarts           = n_restarts,
        n_outer_iters        = n_outer_iters,
    )
    RP = result.expected_obj
    println("  RP = $(round(RP, digits=2))")
end

###############################################################################
# Build the scenario cache and compute the measures
###############################################################################

scenario_cache, rt_mats = build_scenario_cache(data, rt_s, e_s, S)

run_vss_evpi(
    data, rt_s, e_s, S, pi_s, scenario_cache, rt_mats, RP;
    maxTurnaround = maxTurnaround,
    MaxTime_det   = MaxTime_det,
    MaxTime_2nd   = MaxTime_2nd_final,
    top_c         = top_c,
    price_boost   = price_boost,
    hard_penalty  = hard_penalty,
)
