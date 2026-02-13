# Denmark Fictive Data - Vertiports and Population Areas

# Constants
tau = 200
Reasnable_range_factor = 2
Travel_range_to_airort_factor = 1.3
mu = [1.0, 2]
M1 = 1000000
M2 = 1000000

# Vertiport coordinates (15 vertiports across Denmark)
# Roughly mapped to Danish geography (longitude-like, latitude-like)
airport_coords = Dict(
    "V1" => (310, 140),
    "V2" => (220, 125),
    "V3" => (120, 260),
    "V4" => (155, 170),
    "V5" => (100, 100),
    "V6" => (185, 80),
    "V7" => (460, 110),
    "V8" => (80, 155),
    "V9" => (125, 25),
    "V10" => (280, 75)
)


# Operating costs per vertiport
c = Dict(
    "V1" => 12, "V2" => 14, "V3" => 13, "V4" => 25, "V5" => 26,
    "V6" => 15, "V7" => 13, "V8" => 18, "V9" => 19, "V10" => 20
)

# Population areas (8 major population centers)
Population_coords = Dict(
    "P1"  => (470, 105),
    "P2" => (290, 130),
    "P3" => (290, 35),
    "P4" => (90, 120),
    "P5" => (135, 225)
)


# Population demand
Pi = Dict(
    "P1" => 35,    # Copenhagen (largest)
    "P2" => 18,    # Odense (medium)
    "P3" => 22,    # Aarhus (large)
    "P4" => 15,    # Aalborg (medium)
    "P5" => 12,    # Roskilde (medium)
)

# Demand nodes (subset of vertiports for airport connections)
D = ["V1", "V3", "V5"]