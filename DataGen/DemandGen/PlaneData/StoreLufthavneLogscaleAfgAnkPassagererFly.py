# %%
from pathlib import Path
import pandas as pd
import folium
import json
import numpy as np
import branca.colormap as cm

# -------------------------
# 1 Load airport data
# -------------------------

BASE_DIR = Path(__file__).resolve().parent

df = pd.read_excel(BASE_DIR / "EksportSamlet.xlsx", header=[0,1,2])

# Filter only years >= 2020
df = df[df.iloc[:,0] >= 2020]

# Extract terminal departures and arrivals
term_departures = df.xs(('Afg.', 'Term'), level=[1,2], axis=1)
term_arrivals = df.xs(('Ank.', 'Term'), level=[1,2], axis=1)

# Total passengers (Ank + Afg)
term_total = term_departures + term_arrivals
airport_totals = term_total.sum()

# -------------------------
# 2 Airport → municipality
# -------------------------

airport_to_muni = {
    "BILLUND": "Billund Kommune",
    "BORNHOLM/RONNE": "Bornholms Regionskommune",
    "ESBJERG": "Esbjerg Kommune",
    "KOBENHAVN": "Tårnby Kommune",
    "MIDTJYLLANDS LUFTHAVN": "Herning Kommune",
    "ROSKILDE": "Roskilde Kommune",
    "SONDERBORG": "Sønderborg Kommune",
    "THISTED": "Thisted Kommune",
    "AALBORG": "Aalborg Kommune",
    "AARHUS": "Syddjurs Kommune"
}

muni_values = {}

for airport, value in airport_totals.items():

    muni = airport_to_muni.get(airport)

    if muni is not None:
        muni_values[muni] = float(value)

# -------------------------
# 3 Airport coordinates (labels)
# -------------------------

airport_coords = {
    "CPH": (55.6181, 12.6560),
    "BLL": (55.7403, 9.1518),
    "AAL": (57.0928, 9.8492),
    "AAR": (56.3000, 10.6190),
    "RNN": (55.0633, 14.7596),
    "EBJ": (55.5259, 8.5534),
    "RKE": (55.5856, 12.1314),
    "SGD": (54.9644, 9.7917),
    "TED": (57.0688, 8.7052),
    "KRP": (56.3000, 9.2000)
}

# -------------------------
# 4 Load GeoJSON
# -------------------------

with open(BASE_DIR / "kommune.geojson", "r", encoding="utf-8") as f:
    kommuner = json.load(f)

features = [
    f for f in kommuner["features"]
    if f["geometry"]["type"] in ["Polygon", "MultiPolygon"]
]

# -------------------------
# 5 Attach passenger values
# -------------------------

for feature in features:

    name = feature["properties"]["name"]

    value = muni_values.get(name, 0)

    feature["properties"]["passengers"] = f"{int(value):,}"

# -------------------------
# 6 Create map
# -------------------------

m = folium.Map(location=[56.2, 10.0], zoom_start=7)

# -------------------------
# 7 Log-scaled color scale
# -------------------------

log_values = {k: np.log1p(v) for k, v in muni_values.items()}

min_log = min(log_values.values())
max_log = max(log_values.values())

colormap = cm.linear.OrRd_09.scale(min_log, max_log)

# -------------------------
# 8 Style municipalities
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
# 9 Draw municipalities
# -------------------------

folium.GeoJson(
    {"type": "FeatureCollection", "features": features},
    style_function=style_function,
    tooltip=folium.GeoJsonTooltip(
        fields=["name", "passengers"],
        aliases=["Municipality:", "Passengers:"],
        localize=True
    )
).add_to(m)

# -------------------
# -------------------------
# 11 Legend
# -------------------------

colormap.caption = "Total Terminal Passengers (Arrivals + Departures, ≥2020, log scale)"
colormap.add_to(m)

# -------------------------
# 12 Save map
# -------------------------

m.save(BASE_DIR / "airport_municipality_heatmap_post2020.html")

print("Map saved as airport_municipality_heatmap_post2020.html")

# %%



