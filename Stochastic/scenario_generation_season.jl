###############################################################################
# scenario_generation.jl
#
# Computes scenario-dependent travel time and battery consumption parameters
# for the two-stage stochastic eVTOL model.
#
# ─────────────────────────────────────────────────────────────────────────────
# CONFIGURATION — set this before running:
#
const SEASON = :summer   # :summer  or  :winter
#
# :summer  April–September  — more easterly/southerly wind, rarely cold
# :winter  October–March    — SW/W dominant, cold and very cold are realistic
#
# Only the scenario probabilities change between seasons.
# Wind vectors (wx, wy) and temperature multipliers (phi) are identical.
# ─────────────────────────────────────────────────────────────────────────────
#
# Given:
#   - V           : set of vertiport IDs
#   - lat, lon    : coordinate dicts keyed by vertiport ID
#   - rt          : deterministic travel time dict keyed by (i,j)
#   - e           : deterministic battery consumption dict keyed by (i,j)
#   - time_per_km : minutes per km (from params), used to derive v_air
#
# Produces:
#   - rt_s    : Dict{Tuple{Int,Int,Int}, Float64}  keyed by (scenario, i, j)
#   - e_s     : Dict{Tuple{Int,Int,Int}, Float64}  keyed by (scenario, i, j)
#   - S       : 1:length(SCENARIOS)
#   - pi_s    : Dict{Int, Float64}  season-weighted probabilities
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
# Scenario definitions
#
# 13 scenarios derived from DMI wind statistics (2003-2012) for Danish airports.
# Wind vector convention: wx > 0 = eastward, wy > 0 = northward.
# Wind speeds in km/h: 30 ≈ mean from dominant directions, 55 ≈ typical max.
#
# Temperature multiplier phi (battery efficiency degradation):
#   1.00  =>  mild      (> 0 °C)      reference, no penalty
#   1.20  =>  cold      (<= 0 °C)     ~20% extra consumption
#   1.35  =>  very cold (< -8 °C)     ~35% extra consumption
#
# Each scenario has two probability fields:
#   pi_summer : probability April–September (DMI spring/summer wind roses)
#   pi_winter : probability October–March   (DMI autumn/winter wind roses)
#
# Summer: east and south more frequent, strong SW/W less common, cold rare.
# Winter: SW/W strongly dominant, cold and very cold are realistic.
#
# Sources: DMI Technical Report 12-19, wind roses for EKAH/EKBI/EKCH/EKYT.
###############################################################################

const SCENARIOS = [
    # Temperature labels:
    #   AboveZero : > 0°C (covers mild and hot) — phi=1.00, no battery penalty
    #   Cold      : 0°C to -8°C                 — phi=1.20, 20% extra consumption
    #   Very cold : below -8°C                  — phi=1.35, 35% extra consumption
    #
    # "AboveZero" intentionally includes warm/hot summer days: heat-induced
    # battery degradation operates on a multi-day timescale and is not modelled
    # as an operational effect here. phi=1.00 is therefore correct for all
    # above-zero temperatures.
    #
    # Probability breakdown by season:
    #   Summer (Apr-Sep): ~94% AboveZero,  ~6% Cold,  ~0% Very cold
    #   Winter (Oct-Mar): ~55% AboveZero, ~32% Cold, ~13% Very cold
    #
    # Sources: DMI Technical Report 12-19, wind roses EKAH/EKBI/EKCH/EKYT.
    #
    #                            wx      wy    phi   pi_summer  pi_winter
    # ── Calm ─────────────────────────────────────────────────────────────────
    (wx=  0.0, wy=  0.0, phi=1.00, pi_summer=0.16, pi_winter=0.06, label="Calm / AboveZero"),
    (wx=  0.0, wy=  0.0, phi=1.20, pi_summer=0.01, pi_winter=0.03, label="Calm / Cold"),

    # ── Southwest 240° ────────────────────────────────────────────────────────
    # Dominant direction year-round; cold/very cold much more likely in winter.
    # wx=-21, wy=+21 => wind blowing FROM SW TOWARD NE.
    (wx=-21.0, wy= 21.0, phi=1.00, pi_summer=0.16, pi_winter=0.16, label="SW-mod / AboveZero"),
    (wx=-21.0, wy= 21.0, phi=1.20, pi_summer=0.02, pi_winter=0.12, label="SW-mod / Cold"),
    (wx=-21.0, wy= 21.0, phi=1.35, pi_summer=0.00, pi_winter=0.04, label="SW-mod / Very cold"),
    (wx=-39.0, wy= 39.0, phi=1.00, pi_summer=0.08, pi_winter=0.08, label="SW-str / AboveZero"),

    # ── West 270° ─────────────────────────────────────────────────────────────
    # Highest mean wind speed. Strong cold westerlies mainly in winter.
    (wx=-30.0, wy=  0.0, phi=1.00, pi_summer=0.10, pi_winter=0.10, label="W-mod / AboveZero"),
    (wx=-55.0, wy=  0.0, phi=1.20, pi_summer=0.02, pi_winter=0.09, label="W-str / Cold"),
    (wx=-55.0, wy=  0.0, phi=1.35, pi_summer=0.00, pi_winter=0.03, label="W-str / Very cold"),

    # ── East 90° ──────────────────────────────────────────────────────────────
    # Significantly more common in summer (high-pressure over Scandinavia).
    (wx= 30.0, wy=  0.0, phi=1.00, pi_summer=0.24, pi_winter=0.08, label="E-mod / AboveZero"),

    # ── South 180° ────────────────────────────────────────────────────────────
    # More frequent in summer.
    (wx=  0.0, wy= 30.0, phi=1.00, pi_summer=0.20, pi_winter=0.07, label="S-mod / AboveZero"),

    # ── North 360° ────────────────────────────────────────────────────────────
    # Cold northerlies mainly in winter; very cold northerlies are rare.
    (wx=  0.0, wy=-30.0, phi=1.20, pi_summer=0.01, pi_winter=0.08, label="N-mod / Cold"),
    (wx=  0.0, wy=-30.0, phi=1.35, pi_summer=0.00, pi_winter=0.06, label="N-mod / Very cold"),
]

###############################################################################
# Validate both probability sets sum to 1.0 at load time.
###############################################################################

let
    s_sum = sum(sc.pi_summer for sc in SCENARIOS)
    w_sum = sum(sc.pi_winter for sc in SCENARIOS)
    abs(s_sum - 1.0) < 1e-9 || error(
        "Summer probabilities sum to $(s_sum), expected 1.0. " *
        "Please correct the pi_summer fields in the SCENARIOS table."
    )
    abs(w_sum - 1.0) < 1e-9 || error(
        "Winter probabilities sum to $(w_sum), expected 1.0. " *
        "Please correct the pi_winter fields in the SCENARIOS table."
    )
    println("✓ Scenario table loaded: $(length(SCENARIOS)) scenarios")
    println("  pi_summer sums to $(round(s_sum, digits=6)) ✓")
    println("  pi_winter sums to $(round(w_sum, digits=6)) ✓")
end


###############################################################################
# Scenario overview printer
###############################################################################

"""
    print_scenarios()

Print a summary table of all scenarios with both summer and winter
probabilities. Call this at any time to inspect the scenario table
without running the full model.
"""
function print_scenarios()
    println()
    println("=" ^ 75)
    println("SCENARIO TABLE  —  active season: $(SEASON)")
    println("=" ^ 75)
    println(rpad("sc", 4), " | ", rpad("label", 22),
            " | ", lpad("wx", 6), " | ", lpad("wy", 6),
            " | ", lpad("phi", 5),
            " | ", lpad("pi_summer", 9),
            " | ", lpad("pi_winter", 9),
            " | ", lpad("active pi", 9))
    println("-" ^ 75)
    for (sc, scen) in enumerate(SCENARIOS)
        active_pi = SEASON == :summer ? scen.pi_summer : scen.pi_winter
        marker    = active_pi == 0.0 ? "  <- inactive" : ""
        println(rpad(sc, 4), " | ", rpad(scen.label, 22),
                " | ", lpad(scen.wx, 6), " | ", lpad(scen.wy, 6),
                " | ", lpad(scen.phi, 5),
                " | ", lpad(scen.pi_summer, 9),
                " | ", lpad(scen.pi_winter, 9),
                " | ", lpad(round(active_pi, digits=2), 9),
                marker)
    end
    println("-" ^ 75)
    s_sum = sum(sc.pi_summer for sc in SCENARIOS)
    w_sum = sum(sc.pi_winter for sc in SCENARIOS)
    println(rpad("SUM", 4), "   ", rpad("", 22),
            "   ", rpad("", 6), "   ", rpad("", 6), "   ", rpad("", 5),
            " | ", lpad(round(s_sum, digits=4), 9),
            " | ", lpad(round(w_sum, digits=4), 9))
    println("=" ^ 75)
    println()
end

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
    rt_s    : Dict{Tuple{Int,Int,Int}, Float64}  (scenario, i, j) -> travel time
    e_s     : Dict{Tuple{Int,Int,Int}, Float64}  (scenario, i, j) -> battery use
    S       : 1:length(SCENARIOS)
    pi_s    : Dict{Int, Float64}  season-weighted probabilities (from SEASON)
"""
function generate_scenarios(V, lat, lon, rt, e, time_per_km)

    SEASON in (:summer, :winter) || error(
        "SEASON must be :summer or :winter, got :$(SEASON)"
    )

    v_air = 1.0 / time_per_km          # km/min
    v_min = 0.1 * v_air                 # floor to avoid division by zero
    km_per_hour_to_km_per_min = 1.0 / 60.0

    n_scenarios = length(SCENARIOS)
    S = 1:n_scenarios

    # Pick the correct probability column based on SEASON.
    pi_s = Dict(
        sc => (SEASON == :summer ? SCENARIOS[sc].pi_summer : SCENARIOS[sc].pi_winter)
        for sc in S
    )

    rt_s = Dict{Tuple{Int,Int,Int}, Float64}()
    e_s  = Dict{Tuple{Int,Int,Int}, Float64}()

    for (sc, scen) in enumerate(SCENARIOS)
        for i in V, j in V
            i == j && continue

            h = headwind_component(lat[i], lon[i], lat[j], lon[j],
                                   scen.wx, scen.wy)

            # 1 km/h headwind = 1 km/h slower effective airspeed.
            gamma_t = v_air / max(v_min, v_air + h * km_per_hour_to_km_per_min)
            gamma_e = scen.phi * gamma_t

            rt_s[(sc, i, j)] = Int(round(gamma_t * rt[(i, j)]))
            e_s[(sc, i, j)]  = gamma_e * e[(i, j)]
        end
    end

    season_str = SEASON == :summer ? "Summer (Apr–Sep)" : "Winter (Oct–Mar)"
    println("Scenario generation complete:")
    println("  Season         : $(season_str)")
    println("  Scenarios      : $(n_scenarios) × $(length(V) * (length(V) - 1)) route pairs")
    println("  v_air          : $(round(v_air * 60, digits=1)) km/h  ",
            "(v_min = $(round(v_min * 60, digits=1)) km/h)")
    println()
    println(lpad("sc", 4), " | ", rpad("label", 22),
            " | ", lpad("wx", 6), " | ", lpad("wy", 6),
            " | ", lpad("phi", 5), " | ", lpad("pi", 6))
    println("-" ^ 62)
    for (sc, scen) in enumerate(SCENARIOS)
        pi = SEASON == :summer ? scen.pi_summer : scen.pi_winter
        println(lpad(sc, 4), " | ", rpad(scen.label, 22),
                " | ", lpad(scen.wx, 6), " | ", lpad(scen.wy, 6),
                " | ", lpad(scen.phi, 5), " | ", lpad(pi, 6))
    end
    println()

    return rt_s, e_s, S, pi_s
end

# Print scenario table automatically when file is loaded.
print_scenarios()