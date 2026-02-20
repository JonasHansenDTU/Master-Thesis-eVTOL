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
        if name in destination_airports
            if name in charging_airports
                flm.CircleMarker(
                    [lat, lon],
                    radius=8,
                    color="#7a4f01",
                    fill=true,
                    fill_color="#f59e0b",
                    fill_opacity=0.98,
                    popup="Type: Charging + destination airport<br>Name: $(name)"
                ).add_to(m)
            else
                flm.CircleMarker(
                    [lat, lon],
                    radius=7,
                    color="#5b21b6",
                    fill=true,
                    fill_color="#8b5cf6",
                    fill_opacity=0.95,
                    popup="Type: Destination airport<br>Name: $(name)"
                ).add_to(m)
            end
        elseif name in charging_airports
            flm.CircleMarker(
                [lat, lon],
                radius=7,
                color="#14532d",
                fill=true,
                fill_color="#22c55e",
                fill_opacity=0.95,
                popup="Type: Charging airport<br>Name: $(name)"
            ).add_to(m)
        else
            flm.CircleMarker(
                [lat, lon],
                radius=4,
                color="#0f766e",
                fill=true,
                fill_color="#14b8a6",
                fill_opacity=0.92,
                popup="Type: Airport<br>Name: $(name)"
            ).add_to(m)
        end
    end

    for (name, (lat, lon)) in Population_coords
        flm.RegularPolygonMarker(
            [lat, lon],
            number_of_sides=4,
            radius=7,
            rotation=45,
            color="#9a3412",
            fill=true,
            fill_color="#fb923c",
            fill_opacity=0.88,
            popup="Type: Population<br>Name: $(name)"
        ).add_to(m)
    end

    for route in routes
        flm.PolyLine(route, color="#1f77b4", weight=2, opacity=0.8).add_to(m)
    end

    for (pop, airport) in population_links
        pop_lat, pop_lon = Population_coords[pop]
        airport_lat, airport_lon = airport_coords[airport]
        flm.PolyLine(
            [(pop_lat, pop_lon), (airport_lat, airport_lon)],
            color="#2ca02c",
            weight=1,
            opacity=0.6,
            dash_array="5,5"
        ).add_to(m)
    end

    legend_html = """
    <div style="
        position: fixed;
        bottom: 24px;
        left: 24px;
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
        <div><span style="color:#fb923c;">◆</span> Population</div>
        <div><span style="color:#14b8a6;">●</span> Airport</div>
        <div><span style="color:#22c55e;">●</span> Charging airport</div>
        <div><span style="color:#8b5cf6;">●</span> Destination airport</div>
        <div><span style="color:#f59e0b;">●</span> Charging + destination</div>
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
