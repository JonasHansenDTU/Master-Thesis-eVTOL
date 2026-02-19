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
                flm.RegularPolygonMarker(
                    [lat, lon],
                    number_of_sides=4,
                    radius=8,
                    color="#991b1b",
                    fill=true,
                    fill_color="#ef4444",
                    fill_opacity=0.95,
                    popup=name
                ).add_to(m)
            else
                flm.RegularPolygonMarker(
                    [lat, lon],
                    number_of_sides=4,
                    radius=7,
                    color="#1d4ed8",
                    fill=true,
                    fill_color="#3b82f6",
                    fill_opacity=0.9,
                    popup=name
                ).add_to(m)
            end
        elseif name in charging_airports
            flm.CircleMarker(
                [lat, lon],
                radius=6,
                color="#991b1b",
                fill=true,
                fill_color="#ef4444",
                fill_opacity=0.95,
                popup=name
            ).add_to(m)
        else
            flm.CircleMarker(
                [lat, lon],
                radius=4,
                color="#1d4ed8",
                fill=true,
                fill_color="#3b82f6",
                fill_opacity=0.9,
                popup=name
            ).add_to(m)
        end
    end

    for (name, (lat, lon)) in Population_coords
        flm.CircleMarker(
            [lat, lon],
            radius=5,
            color="#c2410c",
            fill=true,
            fill_color="#f97316",
            fill_opacity=0.9,
            popup=name
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

    m.save("map.html")
end


# Data_file_name = "model_data_LonLat.jl" # Change this to "BiggerData.jl" for the larger dataset

# data_file = get(ENV, "MODEL_DATA", joinpath(@__DIR__, "Data", Data_file_name))
# if isfile(data_file)
#     include(data_file)   # defines airport_coords, Population_coords, Pi, tau, etc.
# else
#     error("Data file not found: $data_file")
# end
