using JuMP
using HiGHS

num_city = 3
M = 1000

A = [0 4 2; 4 0 3; 2 3 0]

D = [0 0 0; 2 0 6; 1 7 0]

model = Model(HiGHS.Optimizer)
@variable(model, x[1:num_city, 1:num_city], Bin)
@variable(model, y[1:num_city, 1:num_city]>= 0, Int)
@variable(model, z[1:num_city]>= 0, Int)

@objective(model, Min, sum(x[i,j] for i in 1:num_city for j in 1:num_city))


@constraint(model, [i=1:num_city, j=1:num_city], y[i,j] <= x[i,j]*M)

@constraint(model, [j=1:num_city], sum(y[i,j] for i in 1:num_city) == sum(D[i,j] for i in 1:num_city) + z[j])
@constraint(model, [i=1:num_city], sum(y[i,j] for j in 1:num_city) == sum(D[i,j] for j in 1:num_city) + z[i])


optimize!(model)
if termination_status(model) == OPTIMAL
    println("Optimal solution found!")
    println("Objective value: ", objective_value(model))
    println("x = ", value.(x))
    println("y = ", value.(y))
    println("z = ", value.(z))
else
    println("No optimal solution found. Status: ", termination_status(model))
end