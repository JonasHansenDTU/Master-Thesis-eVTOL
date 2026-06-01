###############################################################################
# scenario_generation.jl
#
# Computes scenario-dependent travel time and battery consumption parameters
# for the two-stage stochastic eVTOL model.
#
# Given:
#   - V         : set of vertiport IDs
#   - lat, lon  : coordinate dicts keyed by vertiport ID
#   - rt        : deterministic travel time dict keyed by (i,j)
#   - e         : deterministic battery consumption dict keyed by (i,j)
#   - time_per_km : minutes per km (from params), used to derive v_air
#
# Produces:
#   - rt_s      : Dict{Tuple{Int,Int,Int}, Float64}  keyed by (scenario, i, j)
#   - e_s       : Dict{Tuple{Int,Int,Int}, Float64}  keyed by (scenario, i, j)
#   - S         : 1:length(SCENARIOS)
#   - pi_s      : Dict{Int, Float64}  empirically weighted probabilities
###############################################################################

###############################################################################
# Geometry helpers
###############################################################################

"""
Compute the forward bearing in degrees from (lat1,lon1) to (lat2,lon2).
Returns a value in [0, 360) where 0 = north, 90 = east.
"""
function bearing_deg(lat1, lon1, lat2, lon2)
    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)
    Δλ = deg2rad(lon2 - lon1)
    x  = sin(Δλ) * cos(φ2)
    y  = cos(φ1) * sin(φ2) - sin(φ1) * cos(φ2) * cos(Δλ)
    θ  = atan(x, y)
    return mod(rad2deg(θ), 360.0)
end

"""
Convert a bearing in degrees to a unit direction vector (dx, dy) where
dx is the eastward component and dy is the northward component.
"""
function bearing_to_unit_vector(b_deg)
    θ  = deg2rad(b_deg)
    dx = sin(θ)   # eastward
    dy = cos(θ)   # northward
    return (dx, dy)
end

"""
Compute the headwind component h^s_ij for the route (i -> j) under wind
vector (wx, wy) in km/h.

Positive h  => headwind (opposing flight, increases time and battery use).
Negative h  => tailwind (assisting flight, decreases time and battery use).

Because d_hat_ij = -d_hat_ji, we correctly get h_ji = -h_ij.
"""
function headwind_component(lat_i, lon_i, lat_j, lon_j, wx, wy)
    b       = bearing_deg(lat_i, lon_i, lat_j, lon_j)
    dx, dy  = bearing_to_unit_vector(b)
    return -(wx * dx + wy * dy)
end

###############################################################################
# Scenario table
#
# 13 scenarios derived from DMI wind statistics (2003-2012) for Danish airports.
# Wind vector convention: wx > 0 = eastward, wy > 0 = northward.
# Wind speeds in km/h: 30 ≈ mean from dominant directions, 55 ≈ typical max.
#
# Temperature multiplier phi (battery efficiency degradation):
#   1.00  =>  mild      (> 0 °C)     reference, no penalty
#   1.20  =>  cold      (<= 0 °C)    ~20% extra consumption
#   1.35  =>  very cold (< -8 °C)    ~35% extra consumption
#
# Probabilities (pi) are empirically weighted from DMI wind roses:
#   SW (240°) is the dominant direction (~15% freq), W (~14%), E+ESE (~15%),
#   S+SSW (~20%), N (~4-5%). Very cold days occur ~5-10 days/year in Denmark.
#   Calm conditions account for ~4-5% of observations.
#   All pi_s sum to 1.0.
###############################################################################

const SCENARIOS = [
    # label                   wx      wy     phi    pi
    # --- calm ---
    (wx=  0.0, wy=  0.0, phi=1.00, pi=0.10, label="Calm / Mild"),
    (wx=  0.0, wy=  0.0, phi=1.20, pi=0.05, label="Calm / Cold"),

    # --- southwest 240° (dominant direction, ~15% freq) ---
    # wx = 30*cos(240°-90°) = 30*sin(240°) ≈ -26, wy = 30*cos(240°) ≈ -15
    # Simplified to equal components for a clean 240° bearing: (-21, 21) means
    # wind blowing FROM southwest TOWARD northeast.
    (wx=-21.0, wy= 21.0, phi=1.00, pi=0.16, label="SW-mod / Mild"),
    (wx=-21.0, wy= 21.0, phi=1.20, pi=0.09, label="SW-mod / Cold"),
    (wx=-21.0, wy= 21.0, phi=1.35, pi=0.04, label="SW-mod / Very cold"),
    (wx=-39.0, wy= 39.0, phi=1.00, pi=0.09, label="SW-str / Mild"),

    # --- west 270° (second most frequent, ~14% freq, highest mean speed) ---
    (wx=-30.0, wy=  0.0, phi=1.00, pi=0.11, label="W-mod / Mild"),
    (wx=-55.0, wy=  0.0, phi=1.20, pi=0.06, label="W-str / Cold"),
    (wx=-55.0, wy=  0.0, phi=1.35, pi=0.03, label="W-str / Very cold"),

    # --- east 90° (E + ESE combined ~15% freq) ---
    (wx= 30.0, wy=  0.0, phi=1.00, pi=0.11, label="E-mod / Mild"),

    # --- south 180° (S + SSW combined ~20% freq) ---
    (wx=  0.0, wy= 30.0, phi=1.00, pi=0.09, label="S-mod / Mild"),

    # --- north 360° (N ~4-5% freq, cold when it blows) ---
    (wx=  0.0, wy=-30.0, phi=1.20, pi=0.05, label="N-mod / Cold"),
    (wx=  0.0, wy=-30.0, phi=1.35, pi=0.02, label="N-mod / Very cold"),
]

# Verify probabilities sum to 1 at load time (fails loudly if table is edited
# incorrectly rather than silently producing wrong results).
const _PI_SUM = sum(sc.pi for sc in SCENARIOS)
abs(_PI_SUM - 1.0) < 1e-9 || error(
    "SCENARIOS probabilities sum to $(_PI_SUM), expected 1.0. " *
    "Please correct the pi fields in the SCENARIOS table."
)

###############################################################################
# Main generation function
###############################################################################

"""
Generate scenario-dependent travel time and battery consumption parameters.

Arguments:
    V           : vector of vertiport IDs
    lat         : Dict{Int, Float64} of latitudes
    lon         : Dict{Int, Float64} of longitudes
    rt          : Dict{Tuple{Int,Int}, Int/Float64} deterministic travel times
    e           : Dict{Tuple{Int,Int}, Float64} deterministic battery consumption
    time_per_km : minutes per km (used to derive nominal airspeed v_air)

Returns:
    rt_s        : Dict{Tuple{Int,Int,Int}, Float64}  (scenario, i, j) -> travel time
    e_s         : Dict{Tuple{Int,Int,Int}, Float64}  (scenario, i, j) -> battery use
    S           : 1:length(SCENARIOS)
    pi_s        : Dict{Int, Float64}  empirically weighted probabilities
"""
function generate_scenarios(V, lat, lon, rt, e, time_per_km)

    v_air = 1.0 / time_per_km          # km/min
    v_min = 0.1 * v_air                 # floor to avoid division by zero
    km_per_hour_to_km_per_min = 1.0 / 60.0

    n_scenarios = length(SCENARIOS)
    S    = 1:n_scenarios

    # Probabilities come from the SCENARIOS table (empirically weighted).
    pi_s = Dict(sc => SCENARIOS[sc].pi for sc in S)

    rt_s = Dict{Tuple{Int,Int,Int}, Float64}()
    e_s  = Dict{Tuple{Int,Int,Int}, Float64}()

    for (sc, scen) in enumerate(SCENARIOS)
        for i in V, j in V
            i == j && continue

            h = headwind_component(lat[i], lon[i], lat[j], lon[j],
                                   scen.wx, scen.wy)

            gamma_t = v_air / max(v_min, v_air + h * km_per_hour_to_km_per_min)
            gamma_e = scen.phi * gamma_t

            rt_s[(sc, i, j)] = Int(round(gamma_t * rt[(i, j)]))
            e_s[(sc, i, j)]  = gamma_e * e[(i, j)]
        end
    end

    println("Scenario generation complete:")
    println("  $(n_scenarios) scenarios × $(length(V) * (length(V) - 1)) route pairs")
    println("  v_air = $(round(v_air * 60, digits=1)) km/h  ",
            "(v_min = $(round(v_min * 60, digits=1)) km/h)")
    println()
    println(lpad("sc", 4), " | ", rpad("label", 20),
            " | ", lpad("wx", 6), " | ", lpad("wy", 6),
            " | ", lpad("phi", 5), " | ", lpad("pi", 6))
    println("-" ^ 60)
    for (sc, scen) in enumerate(SCENARIOS)
        println(lpad(sc, 4), " | ", rpad(scen.label, 20),
                " | ", lpad(scen.wx, 6), " | ", lpad(scen.wy, 6),
                " | ", lpad(scen.phi, 5), " | ", lpad(scen.pi, 6))
    end
    println()

    return rt_s, e_s, S, pi_s
end
