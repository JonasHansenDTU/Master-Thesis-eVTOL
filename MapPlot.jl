using JSON
using PyCall
using GeoInterface
using LibGEOS
using CSV
using DataFrames
using Statistics


function save_folium_map(
    airport_coords,
    Population_coords,
    routes=Vector{Vector{Tuple{Float64,Float64}}}(),
    charging_airports=Vector{String}(),
    population_links=Vector{Tuple{String,String}}(),
    destination_airports=Vector{String}()
)
    flm = pyimport("folium")
    pd = pyimport("pandas")

    latitude_mean = mean([coord[1] for coord in values(airport_coords)])
    longitude_mean = mean([coord[2] for coord in values(airport_coords)])

    m = flm.Map(location=[latitude_mean, longitude_mean], zoom_start=5, tiles="CartoDB positron")

    for (name, (lat, lon)) in airport_coords
        is_destination = name in destination_airports
        is_charging = name in charging_airports
        border_color = is_destination ? "#000000" : (is_charging ? "#7f1d1d" : "#991b1b")
        border_weight = is_destination ? 2 : 1

        if is_charging
            flm.CircleMarker(
                [lat, lon],
                radius=7,
                color=border_color,
                weight=border_weight,
                fill=true,
                fill_color="#ef4444",
                fill_opacity=0.95,
                popup=is_destination ? "Type: Charging airport (Destination)<br>Name: $(name)" : "Type: Charging airport<br>Name: $(name)"
            ).add_to(m)
        else
            flm.CircleMarker(
                [lat, lon],
                radius=4,
                color=border_color,
                weight=border_weight,
                fill=true,
                fill_color="#fca5a5",
                fill_opacity=0.92,
                popup=is_destination ? "Type: Airport (Destination)<br>Name: $(name)" : "Type: Airport<br>Name: $(name)"
            ).add_to(m)
        end

        if is_destination
            flm.Marker(
                [lat, lon],
                icon=flm.DivIcon(
                    html="<div style='width: 14px; height: 14px; display: flex; align-items: center; justify-content: center; font-size: 12px; font-weight: 700; color: #000000;'>✕</div>",
                    icon_size=(14, 14),
                    icon_anchor=(7, 7),
                    class_name="destination-cross"
                )
            ).add_to(m)
        end
    end

    for (name, (lat, lon)) in Population_coords
        flm.RegularPolygonMarker(
            [lat, lon],
            number_of_sides=4,
            radius=5,
            rotation=45,
            color="#7f1d1d",
            fill=true,
            fill_color="#dc2626",
            fill_opacity=0.88,
            popup="Type: Population<br>Name: $(name)"
        ).add_to(m)
    end

    for route in routes
        flm.PolyLine(route, color="#000000", weight=1, opacity=0.70).add_to(m)
    end

    for (pop, airport) in population_links
        pop_lat, pop_lon = Population_coords[pop]
        airport_lat, airport_lon = airport_coords[airport]
        flm.PolyLine(
            [(pop_lat, pop_lon), (airport_lat, airport_lon)],
            color="#000000",
            weight=0.8,
            opacity=0.6,
            dash_array="5,5"
        ).add_to(m)
    end

    legend_html = """
    <div style="
        position: fixed;
        top: 24px;
        right: 24px;
        z-index: 9999;
        background: white;
        border: 1px solid #d1d5db;
        border-radius: 8px;
        padding: 10px 12px;
        font-size: 13px;
        line-height: 1.5;
        box-shadow: 0 2px 8px rgba(0, 0, 0, 0.12);
    ">
        <div style="font-weight: 600; margin-bottom: 6px;">Markers</div>
        <div><span style="color:#dc2626;">◆</span> Population</div>
        <div><span style="color:#fca5a5;">●</span> Airport</div>
        <div><span style="color:#ef4444;">●</span> Charging airport</div>
        <div><span style="color:#111827;">⊗</span> Black border + cross = Destination airport</div>
    </div>
    """
    m.get_root().html.add_child(flm.Element(legend_html))

    m.save("map.html")
end


# Data_file_name = "model_data_LonLat.jl" # Change this to "BiggerData.jl" for the larger dataset

# data_file = get(ENV, "MODEL_DATA", joinpath(@__DIR__, "Data", Data_file_name))
# if isfile(data_file)
#     include(data_file)   # defines airport_coords, Population_coords, Pi, tau, etc.
# else
#     error("Data file not found: $data_file")
# end
