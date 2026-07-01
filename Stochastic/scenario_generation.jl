###############################################################################
# scenario_generation.jl
#
# Fixed set of 13 weather scenarios (wind direction × temperature) for the
# two-stage stochastic eVTOL model. Probabilities are derived as
#
#     pi(scenario) = wind_direction_share × temperature_probability
#
# for the active SEASON, then renormalised over the 13 kept scenarios so they
# sum to 1.
#
# DATA SOURCES
#   • Temperature fractions (AboveZero / Cold / VeryCold) per season are derived
#     from DMI monthly frost-day data 2020–2026 (DMI-DATA.xlsx, Sheet2):
#         frost_pct = Σ frost days / Σ days in season
#         vcold_pct = 25% of frost days in months with absolute min < −8 °C
#         cold_pct  = frost_pct − vcold_pct
#         above_pct = 1 − frost_pct
#   • Wind-direction shares are approximate frequencies read from the DMI
#     wind roses in Technical Report 12-19 (stations EKAH/EKEB/EKOD/EKSB),
#     held constant across seasons. Only the temperature mix changes by season.
#
# SCENARIO REDUCTION
#   The full cross-product is 9 directions × 3 temperatures = 27 combinations.
#   We keep 13 representative ones: all three temperatures for the directions
#   where cold genuinely matters (SW, W, N) and AboveZero only for the
#   mild-dominated directions (Mild, E, S). The probability mass of the omitted
#   combinations is redistributed proportionally via renormalisation. The
#   omitted combinations are individually low-probability and do not materially
#   change the first-stage decision, so 13 trades a small SAA-accuracy loss for
#   a large reduction in second-stage solve cost.
#
# CONFIGURATION:
# SEASON may be overridden by the EVTOL_SEASON environment variable so a driver
# can run different seasons in separate processes without editing this file.
# Defaults to :Sommer when the variable is unset.
const SEASON = Symbol(get(ENV, "EVTOL_SEASON", "Vinter"))    # :Foraar, :Sommer, :Efteraar, :Vinter
#
###############################################################################

###############################################################################
# Geometry helpers
###############################################################################

function bearing_deg(lat1, lon1, lat2, lon2)
    φ1 = deg2rad(lat1); φ2 = deg2rad(lat2); Δλ = deg2rad(lon2 - lon1)
    x = sin(Δλ) * cos(φ2)
    y = cos(φ1) * sin(φ2) - sin(φ1) * cos(φ2) * cos(Δλ)
    return mod(rad2deg(atan(x, y)), 360.0)
end

function bearing_to_unit_vector(b_deg)
    θ = deg2rad(b_deg)
    return (sin(θ), cos(θ))   # (eastward, northward)
end

function headwind_component(lat_i, lon_i, lat_j, lon_j, wx, wy)
    b = bearing_deg(lat_i, lon_i, lat_j, lon_j)
    dx, dy = bearing_to_unit_vector(b)
    return -(wx * dx + wy * dy)
end

###############################################################################
# Temperature mix per season — derived from DMI-DATA.xlsx Sheet2 (2020–2026).
#   (p_aboveZero, p_cold, p_veryCold)   →   phi = 1.00 / 1.20 / 1.35
###############################################################################

const TEMP_PROBS = Dict(
    :Foraar   => (0.831, 0.148, 0.022),  # Marts, April, Maj
    :Sommer   => (1.000, 0.000, 0.000),  # Juni, Juli, August
    :Efteraar => (0.934, 0.055, 0.011),  # September, Oktober, November
    :Vinter   => (0.648, 0.278, 0.074),  # December, Januar, Februar
)

###############################################################################
# Wind-direction shares (approx. from DMI TR12-19 wind roses). Need not sum to
# 1 — normalised inside the probability builder.
###############################################################################

const WIND_SHARES = Dict(
    "N-Mild"        => 0.096,
    "E-Mild"        => 0.14,
    "S-Mild"        => 0.181,
    "W-Mild"        => 0.208,
    "N-Mod"         => 0.028,
    "E-Mod"         => 0.076,
    "S-Mod"         => 0.080,
    "W-Mod"         => 0.162,
    "N-Str"         => 0.002,
    "E-Str"         => 0.004,
    "S-Str"         => 0.004,
    "W-Str"         => 0.019,

)

###############################################################################
# The FIXED 13 scenarios.
#   (wx, wy, phi, wind_key, temp_level, label)
#   temp_level: 1 = AboveZero (phi 1.00), 2 = Cold (phi 1.20), 3 = VeryCold (1.35)
###############################################################################

const SCENARIO_DEFS = [
    #( 0.0,  -9.3, 1.00,  "N-Mild",        1, "N-Mild / AboveZero"),
    #( 0.0,  -9.3, 1.20,  "N-Mild",        2, "N-Mild / Cold"),
    #( 0.0,  -9.3, 1.35,  "N-Mild",        3, "N-Mild / Very Cold"),
    #(-9.3,   0.0, 1.00,  "E-Mild",        1, "E-Mild / AboveZero"),
    #(-9.3,   0.0, 1.20,  "E-Mild",        2, "E-Mild / Cold"),
    #(-9.3,   0.0, 1.35,  "E-Mild",        3, "E-Mild / Very Cold"),
    ( 0.0,   9.3, 1.00,  "S-Mild",        1, "S-Mild / AboveZero"),
    ( 0.0,   9.3, 1.20,  "S-Mild",        2, "S-Mild / Cold"),
    ( 0.0,   9.3, 1.35,  "S-Mild",        3, "S-Mild / Very Cold"),
    ( 9.3,   0.0, 1.00,  "W-Mild",        1, "W-Mild / AboveZero"),
    ( 9.3,   0.0, 1.20,  "W-Mild",        2, "W-Mild / Cold"),
    ( 9.3,   0.0, 1.35,  "W-Mild",        3, "W-Mild / Very Cold"),
    #( 0.0, -27.0, 1.00,  "N-Mod",         1, "N-Mod / AboveZero"),
    #( 0.0, -27.0, 1.20,  "N-Mod",         2, "N-Mod / Cold"),
    #( 0.0, -27.0, 1.35,  "N-Mod",         3, "N-Mod / Very cold"),
    #(-27.0,  0.0, 1.00,  "E-Mod",         1, "E-Mod / AboveZero"),
    #(-27.0,  0.0, 1.20,  "E-Mod",         2, "E-Mod / Cold"),
    #(-27.0,  0.0, 1.35,  "E-Mod",         3, "E-Mod / Very cold"),
    ( 0.0,  27.0, 1.00,  "S-Mod",         1, "S-Mod / AboveZero"),
    ( 0.0,  27.0, 1.20,  "S-Mod",         2, "S-Mod / Cold"),
    #( 0.0,  27.0, 1.35,  "S-Mod",         3, "S-mod / Very cold"),
    ( 27.0,  0.0, 1.00,  "W-Mod",         1, "W-Mod / AboveZero"),
    ( 27.0,  0.0, 1.20,  "W-Mod",         2, "W-Mod / Cold"),
    ( 27.0,  0.0, 1.35,  "W-Mod",         3, "W-mod / Very cold"),
    #( 0.0, -46.0, 1.00,  "N-Str",         1, "N-Str / AboveZero"),
    #( 0.0, -46.0, 1.20,  "N-Str",         2, "N-Str / Cold"),
    #( 0.0, -46.0, 1.35,  "N-Str",         3, "N-Str / Very cold"),
    #(-46.0,  0.0, 1.00,  "E-Str",         1, "E-Str / AboveZero"),
    #(-46.0,  0.0, 1.20,  "E-Str",         2, "E-Str / Cold"),
    #(-46.0,  0.0, 1.35,  "E-Str",         3, "E-Str / Very cold"),
    #( 0.0,  46.0, 1.00,  "S-Str",         1, "S-Str / AboveZero"),
    #( 0.0,  46.0, 1.20,  "S-Str",         2, "S-Str / Cold"),
    #( 0.0,  46.0, 1.35,  "S-Str",         3, "S-Str / Very cold"),
    (46.0,   0.0, 1.00,  "W-Str",         1, "W-Str / AboveZero"),
    #(46.0,   0.0, 1.20,  "W-Str",         2, "W-Str / Cold"),
    #(46.0,   0.0, 1.35,  "W-Str",         3, "W-Str / Very cold"),
]

###############################################################################
# Build the 13 scenarios with season-specific, renormalised probabilities.
###############################################################################

function build_scenarios(season::Symbol)
    season in keys(TEMP_PROBS) || error(
        "SEASON must be :Foraar, :Sommer, :Efteraar or :Vinter, got :$(season)"
    )
    p_temp     = TEMP_PROBS[season]            # (above, cold, vcold)
    total_wind = sum(values(WIND_SHARES))      # normaliser for wind shares

    raw = Float64[]
    for (wx, wy, phi, wkey, tlvl, label) in SCENARIO_DEFS
        push!(raw, (WIND_SHARES[wkey] / total_wind) * p_temp[tlvl])
    end

    total_raw = sum(raw)
    total_raw > 0 || error("All 13 scenarios have zero probability for $(season).")

    scenarios = NamedTuple{(:wx,:wy,:phi,:pi,:label),
                           Tuple{Float64,Float64,Float64,Float64,String}}[]
    for (k, (wx, wy, phi, wkey, tlvl, label)) in enumerate(SCENARIO_DEFS)
        push!(scenarios, (wx=wx, wy=wy, phi=phi,
                          pi = raw[k] / total_raw, label = label))
    end
    return scenarios
end

###############################################################################
# Validation + printer
###############################################################################

function validate_scenarios(scenarios, season)
    total = sum(s.pi for s in scenarios)
    abs(total - 1.0) < 1e-9 || error(
        "$(season) scenario probabilities sum to $(total), expected 1.0"
    )
end

function print_scenarios(scenarios, season)
    println("\n" * "="^72)
    println("SCENARIO TABLE  —  active season: $(season)")
    p_above, p_cold, p_vcold = TEMP_PROBS[season]
    println("  Temperature mix (DMI 2020-2026): AboveZero=$(round(p_above*100,digits=1))%  " *
            "Cold=$(round(p_cold*100,digits=1))%  VCold=$(round(p_vcold*100,digits=1))%")
    println("="^72)
    println(rpad("sc",4), " | ", rpad("label",22), " | ", lpad("wx",6),
            " | ", lpad("wy",6), " | ", lpad("phi",5), " | ", lpad("pi",7))
    println("-"^72)
    for (sc, s) in enumerate(scenarios)
        marker = s.pi == 0.0 ? "  ← zero (skipped)" : (s.pi < 0.01 ? "  ← rare" : "")
        println(rpad(sc,4), " | ", rpad(s.label,22), " | ", lpad(s.wx,6),
                " | ", lpad(s.wy,6), " | ", lpad(s.phi,5),
                " | ", lpad(round(s.pi,digits=4),7), marker)
    end
    println("-"^72)
    println(rpad("SUM",4), "   ", rpad("",22), "   ", rpad("",6), "   ",
            rpad("",6), "   ", rpad("",5), " | ",
            lpad(round(sum(s.pi for s in scenarios),digits=4),7))
    println("="^72 * "\n")
end

###############################################################################
# Main generation function (unchanged interface)
###############################################################################

function generate_scenarios(V, lat, lon, rt, e, time_per_km)
    scenarios = build_scenarios(SEASON)
    validate_scenarios(scenarios, SEASON)

    v_air = 1.0 / time_per_km
    v_min = 0.1 * v_air
    km_per_hour_to_km_per_min = 1.0 / 60.0

    n_scenarios = length(scenarios)
    S    = 1:n_scenarios
    pi_s = Dict(sc => scenarios[sc].pi for sc in S)

    rt_s = Dict{Tuple{Int,Int,Int}, Float64}()
    e_s  = Dict{Tuple{Int,Int,Int}, Float64}()

    for (sc, scen) in enumerate(scenarios)
        for i in V, j in V
            i == j && continue
            h = headwind_component(lat[i], lon[i], lat[j], lon[j], scen.wx, scen.wy)
            gamma_t = v_air / max(v_min, v_air + h * km_per_hour_to_km_per_min)
            gamma_e = scen.phi * gamma_t
            rt_s[(sc, i, j)] = Int(round(gamma_t * rt[(i, j)]))
            e_s[(sc, i, j)]  = gamma_e * e[(i, j)]
        end
    end

    print_scenarios(scenarios, SEASON)
    return rt_s, e_s, S, pi_s
end

# Expose SCENARIOS so downstream code reading SCENARIOS[sc].wx/.wy/.phi/.label
# (the weather-slack helper and the result printers) works unchanged.
const SCENARIOS = build_scenarios(SEASON)

let
    validate_scenarios(SCENARIOS, SEASON)
    print_scenarios(SCENARIOS, SEASON)
end
