function normalize_name(x)
    s = lowercase(strip(String(x)))
    s = replace(s, " " => "_")
    s = replace(s, "-" => "_")
    s = replace(s, "/" => "_")
    s = replace(s, "." => "_")
    s = replace(s, "(" => "")
    s = replace(s, ")" => "")
    return Symbol(s)
end

"Read one Excel sheet into a DataFrame with normalized column names."
function read_sheet(path::String, sheet_name::String)
    df = DataFrame(XLSX.readtable(path, sheet_name, infer_eltypes=true))
    rename!(df, Dict(n => normalize_name(n) for n in names(df)))
    return df
end

"Read first available sheet among candidates into a DataFrame."
function read_sheet_any(path::String, candidates::Vector{String})
    for sheet_name in candidates
        try
            return read_sheet(path, sheet_name)
        catch
            # Try next candidate sheet name.
        end
    end

    error("Could not find any of these sheets: $(candidates)")
end

"Try to find a column among several possible names."
function find_col(df::DataFrame, candidates::Vector{Symbol})
    name_map = Dict(Symbol(String(n)) => n for n in names(df))
    for c in candidates
        if haskey(name_map, c)
            return name_map[c]
        end
    end
    error("Could not find any of these columns in sheet: $(candidates). Found columns: $(names(df))")
end

"Parse coordinate string like '(37.5655, 126.8013)' -> (lat, lon)."
function parse_coordinate_string(s)
    ss = strip(String(s))
    ss = replace(ss, "(" => "")
    ss = replace(ss, ")" => "")
    parts = split(ss, ",")
    if length(parts) != 2
        error("Could not parse coordinate string: $s")
    end
    lat = parse(Float64, strip(parts[1]))
    lon = parse(Float64, strip(parts[2]))
    return lat, lon
end

"Haversine distance in km."
function haversine_km(lat1, lon1, lat2, lon2)
    R = 6371.0
    φ1 = deg2rad(lat1)
    λ1 = deg2rad(lon1)
    φ2 = deg2rad(lat2)
    λ2 = deg2rad(lon2)
    dφ = φ2 - φ1
    dλ = λ2 - λ1
    a = sin(dφ / 2)^2 + cos(φ1) * cos(φ2) * sin(dλ / 2)^2
    c = 2 * atan(sqrt(a), sqrt(1 - a))
    return R * c
end

"""
Build all sets and parameters from Excel + system parameters.
Expected sheets:
    - Infrastructure
    - PassengerGroups 
    - PlaneData
"""
function load_data(excel_file::String, parameter_file::String)

    ###########################################################################
    # Read sheets
    ###########################################################################
    infra = read_sheet(excel_file, "Infrastructure")
    pax   = read_sheet(excel_file, "PassengerGroups")
    plane = read_sheet_any(excel_file, ["PlaneData"])

    ###########################################################################
    # Infrastructure columns
    ###########################################################################
    id_col     = find_col(infra, [:id])
    pads_col   = find_col(infra, [:number_of_parking_pads, :parking_pads, :pads])

    # Coordinates can be either one string column "coordinates"
    # or two numeric columns such as "latitude", "longitude".
    coord_col = if any(Symbol(String(n)) == :coordinates for n in names(infra))
        find_col(infra, [:coordinates])
    else
        nothing
    end
    
    lat_col = if any(Symbol(String(n)) == :latitude for n in names(infra))
        find_col(infra, [:latitude])
    elseif any(Symbol(String(n)) == :lat for n in names(infra))
        find_col(infra, [:lat])
    else
        nothing
    end
    
    lon_col = if any(Symbol(String(n)) == :longitude for n in names(infra))
        find_col(infra, [:longitude])
    elseif any(Symbol(String(n)) == :lon for n in names(infra))
        find_col(infra, [:lon])
    else
        nothing
    end

    ###########################################################################
    # Passenger columns
    ###########################################################################
    group_col  = find_col(pax, [:group, :group_id, :id])
    orig_col   = find_col(pax, [:origin])
    dest_col   = find_col(pax, [:destination])
    time_col   = find_col(pax, [:time, :arrival_time])
    q_col      = find_col(pax, [:number_of_passengers, :passengers, :numberofpassengers])
    stop_col   = find_col(pax, [:stopover_allowed, :stopoverallowed, :stopover])

    ###########################################################################
    # Sets
    ###########################################################################
    V = sort(unique(Int.(infra[!, id_col])))

    plane_id_col = find_col(plane, [:plane_id, :id, :evtol_id, :plane])
    base_vp_col  = find_col(plane, [:base_vertiport, :base_vertiport_id, :base, :home_vertiport])

    N = sort(Int.(plane[!, plane_id_col]))
    if isempty(N)
        error("PlaneData sheet is empty. Add at least one row with Plane ID and Base Vertiport.")
    end
    if length(unique(N)) != length(N)
        error("PlaneData contains duplicate Plane ID values. Each Plane ID must be unique.")
    end

    bv = Dict{Int,Int}()  # base vertiport for each eVTOL
    for r in eachrow(plane)
        n = Int(r[plane_id_col])
        b = Int(r[base_vp_col])
        bv[n] = b
    end

    
    bad_bases = sort([b for b in values(bv) if !(b in V)])
    if !isempty(bad_bases)
        error("PlaneData has invalid Base Vertiport values $(bad_bases). Valid vertiports from Infrastructure are $(V).")
    end

    # Passenger groups
    A = sort(Int.(pax[!, group_col]))

    ###########################################################################
    # System parameters from Table 3
    ###########################################################################
    
    df = XLSX.readtable(parameter_file, "Parameters")
    params = Dict{String, Float64}()

    XLSX.openxlsx(parameter_file) do xf
        sheet = xf["Parameters"]
        
        for row in XLSX.eachrow(sheet)
            key = row[1]
            val = row[2]
            
            if key !== missing && val !== missing
                key_str = String(key)

                # håndter komma som decimalseparator
                val_str = replace(string(val), "," => ".")
                params[key_str] = parse(Float64, val_str)
            end
        end
    end

    fare_direct_per_km    = params["fare_direct_per_km"]
    fare_stopover_factor  = params["fare_stopover_factor"]
    p_penalty             = params["p_penalty"]
    operating_cost_per_km = params["operating_cost_per_km"]
    battery_per_km        = params["battery_per_km"]
    time_per_km           = params["time_per_km"]
    cap_u                 = params["cap_u"]
    opening_cost          = params["opening_cost"]
    bmax                  = params["bmax"]
    bmid                  = params["bmid"]
    b_penalty             = params["b_penalty"]
    bmin                  = params["bmin"]
    ec                    = params["ec"]
    te                    = params["te"]
    w                     = params["w"]
    ET                    = params["ET"]
    M1                    = 1
    M2a                   = 0
    M2b                   = bmax
    M2c                   = bmax + ec * ET
    M3                    = ET

    M = 0:6
    M_no0 = 1:maximum(M)
    M_mid = 1:(maximum(M)-1)
    M_no_last = 0:(maximum(M)-1)

    T = 0:ET
    T_no0 = 1:maximum(T)

    ###########################################################################
    # Node coordinates and parking capacities
    ###########################################################################
    lat = Dict{Int,Float64}()
    lon = Dict{Int,Float64}()
    cap_v = Dict{Int,Int}()

    for r in eachrow(infra)
        j = Int(r[id_col])
        cap_v[j] = Int(r[pads_col])

        if coord_col !== nothing
            lat[j], lon[j] = parse_coordinate_string(r[coord_col])
        else
            if lat_col === nothing || lon_col === nothing
                error("Infrastructure sheet must contain either a 'coordinates' column or latitude/longitude columns.")
            end
            lat[j] = Float64(r[lat_col])
            lon[j] = Float64(r[lon_col])
        end
    end


    ###########################################################################
    # Prices
    ###########################################################################
    prices = read_sheet(parameter_file, "Prices")
    from_col = find_col(prices, [:from])
    to_col   = find_col(prices, [:to])
    fd_sum_col = find_col(prices, [:fd_sum])
    fd_lookup = Dict{Tuple{Int,Int}, Float64}()

    for r in eachrow(prices)
        i = Int(r[from_col])
        j = Int(r[to_col])
        fd_lookup[(i,j)] = Float64(r[fd_sum_col])
    end
    for (i,j) in collect(keys(fd_lookup))
        fd_lookup[(j,i)] = fd_lookup[(i,j)]
    end

    ###########################################################################
    # Derived arc parameters: distance, fd, fs, c, e, rt
    ###########################################################################
    dist = Dict{Tuple{Int,Int},Float64}()
    fd   = Dict{Tuple{Int,Int},Float64}()
    fs   = Dict{Tuple{Int,Int},Float64}()
    c    = Dict{Tuple{Int,Int},Float64}()
    e    = Dict{Tuple{Int,Int},Float64}()
    rt   = Dict{Tuple{Int,Int},Int}()
  

    for i in V, j in V
        dij = haversine_km(lat[i], lon[i], lat[j], lon[j])
        dist[(i,j)] = dij
        fd[(i,j)] = get(fd_lookup, (i,j), 0.0)
        fs[(i,j)] = fare_stopover_factor * fd[(i,j)]
        c[(i,j)]  = operating_cost_per_km * dij
        e[(i,j)]  = battery_per_km * dij
        rt[(i,j)] = Int(ceil(time_per_km * dij))
    end

    ###########################################################################
    # Passenger parameters
    ###########################################################################
    op = Dict{Int,Int}()
    dp = Dict{Int,Int}()
    dt = Dict{Int,Float64}()
    q  = Dict{Int,Int}()
    so = Dict{Int,Int}()
    p  = Dict{Int,Float64}()

    for r in eachrow(pax)
        a = Int(r[group_col])
        op[a] = Int(r[orig_col])
        dp[a] = Int(r[dest_col])
        dt[a] = Float64(r[time_col])
        q[a]  = Int(r[q_col])
        so[a] = Int(r[stop_col])
        p[a]  = p_penalty
    end

    ###########################################################################
    # d[a,i,j] = 1 if passenger group a travels from i to j, else 0
    ###########################################################################
    d = Dict{Tuple{Int,Int,Int},Int}()
    for a in A, i in V, j in V
        d[(a,i,j)] = (i == op[a] && j == dp[a]) ? 1 : 0
    end

    return (
        infra = infra,
        pax = pax,
        plane = plane,
        V = V, A = A, N = collect(N),
        M = collect(M), M_no0 = collect(M_no0), M_mid = collect(M_mid), M_no_last = collect(M_no_last),
        T = collect(T), T_no0 = collect(T_no0),
        bv = bv,
        lat = lat, lon = lon,
        dist = dist, fd = fd_lookup, fs = fs, c = c, e = e, rt = rt,
        op = op, dp = dp, dt = dt, q = q, so = so, p = p, d = d,
        cap_v = cap_v, cap_u = cap_u, opening_cost = opening_cost,
        bmax = bmax, bmid = bmid, b_penalty = b_penalty, bmin = bmin, ec = ec, te = te, w = w, ET = ET, M1 = M1, M2a = M2a, M2b = M2b, M2c = M2c, M3 = M3,
        battery_per_km = battery_per_km
    )
end
