# %%
import pandas as pd
import folium
import json
import numpy as np

# -------------------------
# 1 Load data
# -------------------------

file_path = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"
geo_path = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/kommune.geojson"

df = pd.read_excel(file_path, sheet_name="LTM")
koder = pd.read_excel(file_path, sheet_name="Koder", header=None)
koder.columns = ["kode", "navn"]

with open(geo_path) as f:
    geo = json.load(f)

# -------------------------
# 2 Mapping (kode → Geo navn)
# -------------------------

koder["navn_geo"] = koder["navn"] + " Kommune"
koder.loc[koder["navn"] == "København", "navn_geo"] = "Københavns Kommune"
koder.loc[koder["navn"] == "Bornholm", "navn_geo"] = "Bornholms Regionskommune"

# fjern Christiansø
koder = koder[koder["navn"] != "Christiansø"]

kode_to_geo = dict(zip(koder["kode"], koder["navn_geo"]))

df["FraNavn"] = df["FraKommune"].map(kode_to_geo)
df["TilNavn"] = df["TilKommune"].map(kode_to_geo)

# -------------------------
# 3 GeoJSON → centroid
# -------------------------

features = [
    f for f in geo["features"]
    if f["geometry"]["type"] in ["Polygon", "MultiPolygon"]
]

def get_centroid(feature):
    geom = feature["geometry"]

    if geom["type"] == "Polygon":
        coords = geom["coordinates"][0]
    elif geom["type"] == "MultiPolygon":
        coords = geom["coordinates"][0][0]
    else:
        return None

    lon = sum(p[0] for p in coords) / len(coords)
    lat = sum(p[1] for p in coords) / len(coords)

    return (lat, lon)

kommune_coords = {}

for f in features:
    navn = f["properties"].get("name")
    centroid = get_centroid(f)

    if navn and centroid:
        kommune_coords[navn] = centroid

# -------------------------
# 4 Curved arc funktion
# -------------------------

def bezier_arc(p1, p2, n=60, curvature=0.3):
    lat1, lon1 = p1
    lat2, lon2 = p2

    mid_lat = (lat1 + lat2) / 2
    mid_lon = (lon1 + lon2) / 2

    dlat = lat2 - lat1
    dlon = lon2 - lon1

    ctrl_lat = mid_lat - curvature * dlon
    ctrl_lon = mid_lon + curvature * dlat

    t = np.linspace(0, 1, n)

    lat = (1 - t)**2 * lat1 + 2 * (1 - t) * t * ctrl_lat + t**2 * lat2
    lon = (1 - t)**2 * lon1 + 2 * (1 - t) * t * ctrl_lon + t**2 * lon2

    return list(zip(lat, lon))

# -------------------------
# 5 Filtrér fly flows (ALLE)
# -------------------------

df_fly = df[df["AntalErhvervsture_fly"] > 0].copy()

print(f"Antal fly flows: {len(df_fly)}")

# -------------------------
# 6 Lav kort
# -------------------------

m = folium.Map(location=[56.2, 10.0], zoom_start=7, tiles="CartoDB positron")

# -------------------------
# 7 Kommunegrænser
# -------------------------

folium.GeoJson(
    {"type": "FeatureCollection", "features": features},
    style_function=lambda x: {
        "fillOpacity": 0,
        "color": "#666666",
        "weight": 0.6
    },
    tooltip=folium.GeoJsonTooltip(
        fields=["name"],
        aliases=["Kommune:"]
    )
).add_to(m)

# -------------------------
# 8 Tegn fly flows (curved)
# -------------------------

max_flow = df_fly["AntalErhvervsture_fly"].max()

for _, row in df_fly.iterrows():

    fra = row["FraNavn"]
    til = row["TilNavn"]
    value = row["AntalErhvervsture_fly"]

    if fra in kommune_coords and til in kommune_coords:

        arc = bezier_arc(
            kommune_coords[fra],
            kommune_coords[til],
            curvature=0.25
        )

        weight = 1 + (value / max_flow) * 6

        folium.PolyLine(
            arc,
            weight=weight,
            color="#6A5ACD",  # lilla/blå
            opacity=0.7,
            tooltip=f"{fra} → {til}<br>{value:.1f} flyture"
        ).add_to(m)

# -------------------------
# 9 Marker noder (lufthavne-ish kommuner)
# -------------------------

active_nodes = set(df_fly["FraNavn"]).union(df_fly["TilNavn"])

for navn in active_nodes:
    if navn in kommune_coords:
        folium.CircleMarker(
            location=kommune_coords[navn],
            radius=5,
            color="black",
            fill=True,
            fill_color="#6A5ACD",
            fill_opacity=0.9,
            tooltip=navn
        ).add_to(m)

# -------------------------
# 10 Gem kort
# -------------------------

output_file = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/flow_map_fly_curved.html"

m.save(output_file)

print(f"Kort gemt: {output_file}")