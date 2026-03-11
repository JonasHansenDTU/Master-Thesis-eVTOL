# %%
from pathlib import Path

import pandas as pd
import folium
import json
import numpy as np
import branca.colormap as cm

# -------------------------
# 1 Load passenger data
# -------------------------

BASE_DIR = Path(__file__).resolve().parent

df = pd.read_excel(BASE_DIR / "Passagerer_clean.xlsx")

# Filter only years >= 2020
df = df[df["År"] >= 2020]

# Keep only Indenrigs Total rows
df = df[
    (df["indenrigs Udenrigs"] == "Indenrigs") &
    (df["Flyvning"] == "Total")
]

# Airport columns start after the first 3 columns
airport_data = df.iloc[:,3:]

# Sum passengers per airport
airport_totals = airport_data.sum()

# -------------------------
# 2 Airport → municipality
# -------------------------

airport_to_muni = {
    "AERO": "Ærø Kommune",
    "HERNING": "Herning Kommune",
    "KOBENHAVN VANDFLYVEPLADS (WATER AD)": "Københavns Kommune",
    "KOLDING/VAMDRUP": "Kolding Kommune",
    "LOLLAND FALSTER/MARIBO": "Lolland Kommune",
    "ODENSE/HANS CHRISTIAN ANDERSEN AIRPORT": "Nordfyns Kommune",
    "RANDERS": "Randers Kommune",
    "SINDAL": "Hjørring Kommune",
    "STAUNING": "Ringkøbing-Skjern Kommune",
    "TONDER": "Tønder Kommune",
    "TAASINGE/ELVIRA MADIGAN AIRPORT": "Svendborg Kommune",
    "VIBORG": "Viborg Kommune",
    "VOJENS/SKRYDSTRUP (MIL)": "Haderslev Kommune",
    "AARHUS VANDFLYVEPLADS (WATER AD)": "Aarhus Kommune"
}

muni_values = {}
muni_airport = {}

for airport, value in airport_totals.items():

    muni = airport_to_muni.get(airport)

    if muni is not None:
        muni_values[muni] = float(value)
        muni_airport[muni] = airport

# -------------------------
# 3 Load GeoJSON
# -------------------------

with open(BASE_DIR / "kommune.geojson") as f:
    kommuner = json.load(f)

features = [
    f for f in kommuner["features"]
    if f["geometry"]["type"] in ["Polygon", "MultiPolygon"]
]

# -------------------------
# 4 Attach airport + passengers
# -------------------------

for feature in features:

    name = feature["properties"]["name"]

    value = muni_values.get(name, 0)
    airport = muni_airport.get(name, "None")

    feature["properties"]["airport"] = airport
    feature["properties"]["passengers"] = f"{int(value):,}"

# -------------------------
# 5 Create map
# -------------------------

m = folium.Map(location=[56.2, 10.0], zoom_start=7)

# -------------------------
# 6 Log-scaled color scale
# -------------------------

log_values = {k: np.log1p(v) for k, v in muni_values.items()}

min_log = min(log_values.values())
max_log = max(log_values.values())

colormap = cm.linear.OrRd_09.scale(min_log, max_log)

# -------------------------
# 7 Style municipalities
# -------------------------

def style_function(feature):

    name = feature["properties"]["name"]

    value = muni_values.get(name)

    if value is None:
        return {
            "fillOpacity": 0,
            "color": "black",
            "weight": 1
        }

    return {
        "fillColor": colormap(np.log1p(value)),
        "color": "black",
        "weight": 1,
        "fillOpacity": 0.9
    }

# -------------------------
# 8 Draw municipalities
# -------------------------

folium.GeoJson(
    {"type": "FeatureCollection", "features": features},
    style_function=style_function,
    tooltip=folium.GeoJsonTooltip(
        fields=["name", "airport", "passengers"],
        aliases=["Municipality:", "Airport:", "Passengers:"],
        localize=True
    )
).add_to(m)

# -------------------------
# 9 Legend
# -------------------------

colormap.caption = "Domestic passengers ≥2020 (log scale)"
colormap.add_to(m)

# -------------------------
# 10 Save map
# -------------------------

m.save(BASE_DIR / "domestic_airfields_heatmap.html")

print("Map saved as domestic_airfields_heatmap.html")


