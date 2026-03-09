# ...existing code...
# Define only data (no logic). Example:
tau = 300
Reasnable_range_factor = 2
Travel_range_to_airort_factor = 1.74
mu = [1.0, 0.5]
M1 = 1000
M2 = 1000

airport_coords = Dict(
    "A" => (55.6761, 12.5683),
    "B" => (56.1629, 10.2039),
    "C" => (55.3959, 10.4024),
    "D" => (55.7400, 9.1510)
)

c = Dict("A"=>10, "B"=>20, "C"=>15, "D"=>25)

Population_coords = Dict(
    "K1"=>(55.5761, 12.6683), 
    "K2"=>(56.2629, 10.1039)
)

Pi = Dict("K1"=>1000, "K2"=>1500)

D = ["D"]
# any other constants / lists you need (E, P start empty in main script)
# ...existing code...