using XLSX
using HTTP
using JSON3

file = "inputData/Parameters.xlsx"
sheet_name = "Prices"

# -----------------------------
# Parsing helpers
# -----------------------------

function parse_number(x)
    if ismissing(x) || x === nothing || strip(string(x)) == ""
        return nothing
    end

    return parse(Float64, replace(strip(string(x)), "," => "."))
end

function parse_coord(s)
    if ismissing(s) || s === nothing || strip(string(s)) == ""
        return nothing
    end

    txt = replace(string(s), "(" => "", ")" => "")
    parts = split(txt, ",")

    if length(parts) != 2
        error("Could not parse coordinate: $s")
    end

    lat = parse_number(parts[1])
    lon = parse_number(parts[2])

    return lat, lon
end

# -----------------------------
# Distance calculation
# -----------------------------

function haversine_km(lat1, lon1, lat2, lon2)
    r = 6371.0

    φ1 = deg2rad(lat1)
    φ2 = deg2rad(lat2)
    Δφ = deg2rad(lat2 - lat1)
    Δλ = deg2rad(lon2 - lon1)

    a = sin(Δφ / 2)^2 + cos(φ1) * cos(φ2) * sin(Δλ / 2)^2
    c = 2 * asin(sqrt(a))

    return r * c
end

# -----------------------------
# OSRM route data
# -----------------------------

function osrm_route(lat1, lon1, lat2, lon2)
    # OSRM bruger lon,lat — ikke lat,lon
    url = "http://router.project-osrm.org/route/v1/driving/$(lon1),$(lat1);$(lon2),$(lat2)?overview=false"

    resp = HTTP.get(url)
    data = JSON3.read(String(resp.body))

    if data.code != "Ok"
        error("OSRM error: $(data.code)")
    end

    duration_minutes = data.routes[1].duration / 60.0
    road_distance_km = data.routes[1].distance / 1000.0

    return duration_minutes, road_distance_km
end

function route_key(lat1, lon1, lat2, lon2)
    return (
        round(lat1, digits=6),
        round(lon1, digits=6),
        round(lat2, digits=6),
        round(lon2, digits=6)
    )
end

function route_minutes_and_km(lat1, lon1, lat2, lon2, cache)
    key = route_key(lat1, lon1, lat2, lon2)

    if haskey(cache, key)
        return cache[key]
    end

    result = osrm_route(lat1, lon1, lat2, lon2)
    cache[key] = result

    # Vær sød mod den offentlige OSRM-server
    sleep(0.2)

    return result
end

# -----------------------------
# Excel helpers
# -----------------------------

function find_header(headers, name)
    idx = findfirst(x -> strip(string(x)) == name, headers)

    if idx === nothing
        error("Could not find column: $name")
    end

    return idx
end

function find_optional_header(headers, name)
    return findfirst(x -> strip(string(x)) == name, headers)
end

function read_headers(sheet, last_col)
    return [sheet[1, c] for c in 1:last_col]
end

function delete_column_if_exists!(sheet, headers, last_row, last_col, column_name)
    col_idx = find_optional_header(headers, column_name)

    if col_idx === nothing
        return last_col
    end

    for r in 1:last_row
        for c in col_idx:(last_col - 1)
            sheet[r, c] = sheet[r, c + 1]
        end
        sheet[r, last_col] = ""
    end

    return last_col - 1
end

function ensure_column_after!(sheet, headers, last_row, last_col, new_column_name, after_column_name)
    existing_idx = find_optional_header(headers, new_column_name)

    if existing_idx !== nothing
        return existing_idx, last_col
    end

    after_idx = find_header(headers, after_column_name)
    insert_idx = after_idx + 1

    for r in 1:last_row
        for c in last_col:-1:insert_idx
            sheet[r, c + 1] = sheet[r, c]
        end
        sheet[r, insert_idx] = ""
    end

    sheet[1, insert_idx] = new_column_name

    return insert_idx, last_col + 1
end

function find_parameter_value(sheet, parameter_name)
    dim = XLSX.get_dimension(sheet)

    for r in dim.start.row_number:dim.stop.row_number
        for c in dim.start.column_number:dim.stop.column_number
            val = sheet[r, c]

            if !ismissing(val) && val !== nothing && strip(string(val)) == parameter_name
                parameter_value = sheet[r, c + 1]
                parsed = parse_number(parameter_value)

                if parsed === nothing
                    error("Could not parse value for parameter: $parameter_name")
                end

                return parsed
            end
        end
    end

    error("Could not find parameter: $parameter_name")
end

# -----------------------------
# Main script
# -----------------------------

cache = Dict()

XLSX.openxlsx(file, mode="rw") do xf
    sheet = xf[sheet_name]

    dim = XLSX.get_dimension(sheet)
    last_row = dim.stop.row_number
    last_col = dim.stop.column_number

    headers = read_headers(sheet, last_col)

    # Delete the old unused "fd max" column if it exists.
    last_col = delete_column_if_exists!(sheet, headers, last_row, last_col, "fd max")
    headers = read_headers(sheet, last_col)

    # Add road-distance column between "cord_to" and "Bil tid (min)" if it does not already exist.
    _, last_col = ensure_column_after!(sheet, headers, last_row, last_col, "Bil km", "cord_to")
    headers = read_headers(sheet, last_col)

    col_from        = find_header(headers, "From")
    col_to          = find_header(headers, "To")
    col_cord_from   = find_header(headers, "cord_from")
    col_cord_to     = find_header(headers, "cord_to")
    col_road_km     = find_header(headers, "Bil km")
    col_car_time    = find_header(headers, "Bil tid (min)")
    col_km          = find_header(headers, "Km fugleflugt")
    col_evtol_time  = find_header(headers, "eVTOL flyvetid (min)")
    col_faster      = find_header(headers, "Hvor meget hurtigere")
    col_fd_km       = find_header(headers, "Fd km")
    col_fd_time     = find_header(headers, "Fd tid")
    col_fd_sum      = find_header(headers, "fd_sum")

    time_per_km        = find_parameter_value(sheet, "time_per_km")
    passenger_km_price = find_parameter_value(sheet, "Passager km pris")
    time_value         = find_parameter_value(sheet, "Tidsværdi")
    lambda             = find_parameter_value(sheet, "Lambda")
    base_fee           = find_parameter_value(sheet, "Base fee")

    println("Parameters:")
    println("time_per_km = $time_per_km")
    println("passenger_km_price = $passenger_km_price")
    println("time_value = $time_value")
    println("lambda = $lambda")
    println("base_fee = $base_fee")

    for r in 2:last_row
        from_val = sheet[r, col_from]
        to_val   = sheet[r, col_to]

        if ismissing(from_val) || ismissing(to_val) || from_val === nothing || to_val === nothing
            continue
        end

        if from_val == to_val
            sheet[r, col_road_km] = 0
            sheet[r, col_car_time] = 0
            sheet[r, col_km] = 0
            sheet[r, col_evtol_time] = 0
            sheet[r, col_faster] = 0
            sheet[r, col_fd_km] = 0
            sheet[r, col_fd_time] = 0
            sheet[r, col_fd_sum] = 0
            continue
        end

        coord_from = parse_coord(sheet[r, col_cord_from])
        coord_to   = parse_coord(sheet[r, col_cord_to])

        if coord_from === nothing || coord_to === nothing
            continue
        end

        lat1, lon1 = coord_from
        lat2, lon2 = coord_to

        try
            route_minutes, route_km = route_minutes_and_km(lat1, lon1, lat2, lon2, cache)

            car_time = round(route_minutes, digits=0)
            road_km = round(route_km, digits=0)
            km = round(haversine_km(lat1, lon1, lat2, lon2), digits=0)

            evtol_time = round(km * time_per_km, digits=1)
            faster = evtol_time == 0 ? 0 : round(car_time / evtol_time, digits=1)

            fd_km = round(km * passenger_km_price + base_fee, digits=0)
            fd_time = round((car_time - evtol_time) * time_value * lambda, digits=0)
            fd_sum = round(fd_km + fd_time, digits=0)

            sheet[r, col_road_km] = road_km
            sheet[r, col_car_time] = car_time
            sheet[r, col_km] = km
            sheet[r, col_evtol_time] = evtol_time
            sheet[r, col_faster] = faster
            sheet[r, col_fd_km] = fd_km
            sheet[r, col_fd_time] = fd_time
            sheet[r, col_fd_sum] = fd_sum

            println(
                "Row $r: $from_val → $to_val | " *
                "road_km=$(road_km), car=$(car_time) min, fugleflugt_km=$(km), evtol=$(evtol_time), faster=$(faster), fd_sum=$(fd_sum)"
            )

        catch e
            println("Could not calculate row $r: $from_val → $to_val")
            println(e)
        end
    end
end
