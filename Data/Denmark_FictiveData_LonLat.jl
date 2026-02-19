# Denmark Fictive Data - Vertiports and Population Areas

# Constants
tau = 250
Reasnable_range_factor = 2
Travel_range_to_airort_factor = 1.3
mu = [1.0, 1.5]
M1 = 1000000
M2 = 1000000

# Vertiport coordinates (lat, lon) for Danish cities
airport_coords = Dict(
    "V1" => (55.67519054505104, 12.5381176956762),  # Copenhagen
    "V2" => (56.16246384036334, 10.170032042767225),  # Aarhus
    "V3" => (55.39648487553288, 10.400276215511703),  # Odense
    "V4" => (57.027254122761725, 9.9499495483339),   # Aalborg
    "V5" => (55.487271260981885, 8.495563227580453),   # Esbjerg
    "V6" => (56.45816752400199, 10.065278419943223),  # Randers
    "V7" => (55.4957202757602, 9.46278579491893),   # Kolding
    "V8" => (55.86134265390743, 9.847812369064993),   # Horsens
    "V9" => (55.71362481895501, 9.5207891555811),   # Vejle
    "V10" => (55.64199391085855, 12.069378458129927),  # Roskilde
    "V11" => (54.91118798381618, 9.81451756940919)  # Sønderborg
)


# Operating costs per vertiport
c = Dict(
    "V1" => 12, "V2" => 14, "V3" => 13, "V4" => 25, "V5" => 26,
    "V6" => 15, "V7" => 13, "V8" => 18, "V9" => 19, "V10" => 20, "V11" => 22
)

# Population areas (lat, lon) for major population centers
Population_coords = Dict(
    "P1" => (55.85399258428812, 12.381507893944194),  # Nordsjældand
    "P2" => (55.32797822347244, 10.536651989340676),  # Fyn
    "P3" => (57.22854526249185, 9.863024883325647),  # Nordjylland
    "P4" => (55.02534591045084, 9.371776929993308),   # Aabenraa
    "P5" => (56.154939809259844, 9.55265617771332),   # Silkeborg
    "P6" => (55.200175371546216, 11.847819233431546)   # Næstved
)


# Population demand
Pi = Dict(
    "P1" => 35,    # Nordsjældand (largest)
    "P2" => 18,    # Fyn (medium)
    "P3" => 22,    # Nordjylland (large)
    "P4" => 15,    # Sønderborg (medium)
    "P5" => 12,    # Silkeborg (medium)
    "P6" => 10,    # Næstved (small)
)

# Demand nodes (subset of vertiports for airport connections)
D = ["V1", "V3", "V5"]