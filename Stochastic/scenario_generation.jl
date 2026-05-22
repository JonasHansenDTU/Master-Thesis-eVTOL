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
#   - scenarios : Vector of NamedTuples with fields (wx, wy, phi, label)
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
    θ  = atan(x, y)   # note: atan(y,x) in most languages; Julia atan(x,y) = atan2
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

Because d_hat_ij = -d_hat_ji, we correctly get h_ji = -h_ij, capturing
the physical asymmetry between outbound and return legs.
"""
function headwind_component(lat_i, lon_i, lat_j, lon_j, wx, wy)
    b       = bearing_deg(lat_i, lon_i, lat_j, lon_j)
    dx, dy  = bearing_to_unit_vector(b)
    # h = -(w · d_hat) = -(wx*dx + wy*dy)
    return -(wx * dx + wy * dy)
end

###############################################################################
# Scenario table
###############################################################################
# 18 scenarios: 9 wind conditions × 2 temperature levels (mild, cold).
#
# Wind vector convention:
#   wx > 0  =>  wind blowing eastward
#   wy > 0  =>  wind blowing northward
#
# Wind speeds: 0, 30 (moderate), 60 (strong) km/h.
#
# Temperature multiplier phi:
#   1.00  =>  mild  (> 0 °C), no battery penalty
#   1.20  =>  cold  (~-10 °C), 20% extra battery consumption
###############################################################################

const SCENARIOS = [
    # Wind still
    (wx=  0.0, wy=  0.0, phi=1.00, label="Still / Mild"),
    (wx=  0.0, wy=  0.0, phi=1.00, label="Still / Mild"),
    (wx=  0.0, wy=  0.0, phi=1.20, label="Still / Cold"),
    # Moderate north wind  (blowing northward, wy > 0)
    (wx=  0.0, wy= 30.0, phi=1.00, label="N-30 / Mild"),
    (wx=  0.0, wy= 30.0, phi=1.20, label="N-30 / Cold"),
    # Moderate south wind  (blowing southward, wy < 0)
    (wx=  0.0, wy=-30.0, phi=1.00, label="S-30 / Mild"),
    (wx=  0.0, wy=-30.0, phi=1.20, label="S-30 / Cold"),
    # # Moderate north wind  (blowing northward, wy > 0)
    (wx=  0.0, wy= 30.0, phi=1.00, label="N-30 / Mild"),
    (wx=  0.0, wy= 30.0, phi=1.20, label="N-30 / Cold"),
    # # Moderate south wind  (blowing southward, wy < 0)
    # (wx=  0.0, wy=-30.0, phi=1.00, label="S-30 / Mild"),
    # (wx=  0.0, wy=-30.0, phi=1.20, label="S-30 / Cold"),
    # # Moderate east wind   (blowing eastward,  wx > 0)
    (wx= 30.0, wy=  0.0, phi=1.00, label="E-30 / Mild"),
    (wx= 30.0, wy=  0.0, phi=1.20, label="E-30 / Cold"),
    # # Moderate west wind   (blowing westward,  wx < 0)
    (wx=-30.0, wy=  0.0, phi=1.00, label="W-30 / Mild"),
    (wx=-30.0, wy=  0.0, phi=1.20, label="W-30 / Cold"),
    # # Strong north wind
    (wx=  0.0, wy= 60.0, phi=1.00, label="N-60 / Mild"),
    (wx=  0.0, wy= 60.0, phi=1.20, label="N-60 / Cold"),
    # # Strong south wind
    (wx=  0.0, wy=-60.0, phi=1.00, label="S-60 / Mild"),
    (wx=  0.0, wy=-60.0, phi=1.20, label="S-60 / Cold"),
    # # Strong east wind
    (wx= 60.0, wy=  0.0, phi=1.00, label="E-60 / Mild"),
    (wx= 60.0, wy=  0.0, phi=1.20, label="E-60 / Cold"),
    # # Strong west wind
    (wx=-60.0, wy=  0.0, phi=1.00, label="W-60 / Mild"),
    (wx=-60.0, wy=  0.0, phi=1.20, label="W-60 / Cold"),
]

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
    S           : 1:length(SCENARIOS)  (scenario index set)
    pi_s        : Dict{Int, Float64}   equal probabilities 1/|S|
"""
function generate_scenarios(V, lat, lon, rt, e, time_per_km)

    # Nominal airspeed in km/min; minimum set to 10% of nominal to avoid
    # division by zero or physically impossible values under extreme headwinds.
    v_air = 1.0 / time_per_km
    v_min = 0.1 * v_air
    km_per_hour_to_km_per_min = 1.0 / 60.0

    n_scenarios = length(SCENARIOS)
    S    = 1:n_scenarios
    pi_s = Dict(sc => 1.0 / n_scenarios for sc in S)

    rt_s = Dict{Tuple{Int,Int,Int}, Float64}()
    e_s  = Dict{Tuple{Int,Int,Int}, Float64}()

    for (sc, scen) in enumerate(SCENARIOS)
        for i in V, j in V
            

            h = headwind_component(lat[i], lon[i], lat[j], lon[j],
                                   scen.wx, scen.wy)

            # Travel-time multiplier:
            #   gamma_t = v_air / max(v_min, v_air + h/60)
            # h > 0 (headwind)  => effective airspeed decreases => gamma_t > 1
            # h < 0 (tailwind)  => effective airspeed increases => gamma_t < 1
            gamma_t = v_air / max(v_min, v_air + h * km_per_hour_to_km_per_min)

            # Battery multiplier combines wind effect (same aerodynamic origin
            # as travel time) with temperature effect (phi, route-independent).
            gamma_e = scen.phi * gamma_t

            

            rt_s[(sc, i, j)] = Int(round(0.9 * rt[(i, j)]))
            rt_s[(sc, i, j)] = Int(round(gamma_t * rt[(i, j)]))
            e_s[(sc, i, j)]  = gamma_e * e[(i, j)]
        end
    end

    println("Scenario generation complete:")
    println("  $(n_scenarios) scenarios × $(length(V)*(length(V)-1)) route pairs")
    println("  v_air = $(round(v_air * 60, digits=1)) km/h  ",
            "(v_min = $(round(v_min * 60, digits=1)) km/h)")
    println()
    println(lpad("sc", 4), " | ", lpad("label", 14),
            " | ", lpad("wx", 6), " | ", lpad("wy", 6),
            " | ", lpad("phi", 5))
    println("-" ^ 48)
    for (sc, scen) in enumerate(SCENARIOS)
        println(lpad(sc, 4), " | ", rpad(scen.label, 14),
                " | ", lpad(scen.wx, 6), " | ", lpad(scen.wy, 6),
                " | ", lpad(scen.phi, 5))
    end
    println()

    return rt_s, e_s, S, pi_s
end
