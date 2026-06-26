include("MainCall_stochastic.jl")  # or however you load data + scenarios

# After data, rt_s, e_s, S, pi_s are loaded, build a first-stage decision with
# 1 active eVTOL and NO committed passengers (the case that fails):
test_fsd = make_fsd(Set{Int}(), Set{Int}([1]), data,
                    allPlaneSolution([planeSolution(Int32(0), Int32[Int32(data.bv[n])], Int32[]) for n in data.N]))

# Build scenario 1's data and solve its second stage directly:
demand_realized = make_demand_realizations(data, S; rho=0.5, seed=20260101)
cache, rtmats = build_scenario_cache(data, rt_s, e_s, S; demand_realized=demand_realized)
sc = first(S)
r = solve_second_stage(test_fsd, cache[sc], rtmats[sc];
                       maxTurnaround=100, MaxTime_2nd=Int32(15),
                       top_c=4, price_boost=8.0, hard_penalty=50000.0)
println("Pool B available in sc$sc: ", cache[sc].B_realized)
println("Served: ", [ass.group for ass in r.assignments])
println("Profit: ", r.profit)
print_chromosome_table(r.sol)