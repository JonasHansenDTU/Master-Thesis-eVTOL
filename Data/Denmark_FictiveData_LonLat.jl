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
    "V11" => (54.91118798381618, 9.81451756940919),  # Sønderborg
    # "V12" => (56.036004, 12.613525),   # Helsingør
    "V13" => (55.403756, 11.354718),   # Slagelse
    "V14" => (55.059819, 10.606294),   # Svendborg
    # "V15" => (56.494720, 8.976280),    # Holstebro
    # "V16" => (56.565987, 9.028457),    # Herning
    # "V17" => (56.706861, 11.539428),   # Hjørring
    # "V18" => (54.769092, 11.874978),   # Nykøbing Falster
    # "V19" => (55.770485, 12.503779),   # Hillerød
     "V20" => (55.0985561, 14.699377),    # Rønne
    # "V21" => (55.709270, 11.722110)    # Kalundborg
)


# Operating costs per vertiport
c = Dict(
    "V1" => 12, "V2" => 14, "V3" => 13, "V4" => 25, "V5" => 26,
    "V6" => 15, "V7" => 13, "V8" => 18, "V9" => 19, "V10" => 20, "V11" => 22,
    # "V12" => 16,
    "V13" => 17,
    "V14" => 15,
    # "V15" => 18,
    # "V16" => 19,
    # "V17" => 20,
    # "V18" => 21,
    # "V19" => 16,
     "V20" => 14,
    # "V21" => 18
)

# Population areas (lat, lon) for major population centers
Population_coords = Dict(
    "P1" => (55.85399258428812, 12.381507893944194),  # Nordsjældand
    "P2" => (55.32797822347244, 10.536651989340676),  # Fyn
    "P3" => (57.22854526249185, 9.863024883325647),  # Nordjylland
    "P4" => (55.02534591045084, 9.371776929993308),   # Aabenraa
    "P5" => (56.154939809259844, 9.55265617771332),   # Silkeborg
    "P6" => (55.200175371546216, 11.847819233431546),   # Næstved
    "P7" => (56.162939, 8.977064),    # Midtjylland West
    "P8" => (55.350429, 11.167570),   # Korsør
    "P9" => (54.938724, 9.919396),    # Haderslev
    "P10" => (56.795236, 8.863991),   # Thisted
    "P11" => (55.100152, 14.706270),  # Bornholm
    "P12" => (55.942047, 9.127533),   # Viborg
    # "P13" => (55.490400, 11.343798),  # Ringsted
    # "P14" => (56.360091, 10.033400),  # Skanderborg
    # "P15" => (55.251799, 12.301248),  # Køge
    # "P16" => (56.091006, 12.022500)   # Frederikssund
)


# Population demand
Pi = Dict(
    "P1" => 35,    # Nordsjældand (largest)
    "P2" => 18,    # Fyn (medium)
    "P3" => 22,    # Nordjylland (large)
    "P4" => 15,    # Sønderborg (medium)
    "P5" => 12,    # Silkeborg (medium)
    "P6" => 10,    # Næstved (small)
    "P7" => 14,
    "P8" => 9,
    "P9" => 11,
    "P10" => 8,
    "P11" => 7,
    "P12" => 13,
    # "P13" => 10,
    # "P14" => 12,
    # "P15" => 11,
    # "P16" => 9
)

# Demand nodes (subset of vertiports for airport connections)
D = ["V1", "V3", "V5", "V4"]