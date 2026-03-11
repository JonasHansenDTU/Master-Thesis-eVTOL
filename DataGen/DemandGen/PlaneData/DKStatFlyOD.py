# %% [markdown]
# Only 2025 data
# 

# %%
from unittest.mock import Base

from matplotlib.path import Path
import pandas as pd
import folium
import numpy as np
import branca.colormap as cm
from pathlib import Path

# -------------------------
# 1 Load OD matrix
# -------------------------

BASE_DIR = Path(__file__).resolve().parent

df = pd.read_excel(BASE_DIR / "DKStatFly.xlsx", index_col=0)
df = df.apply(pd.to_numeric, errors="coerce").fillna(0)

# -------------------------
# 2 Airport coordinates
# -------------------------

airport_coords = {
    "København": (55.6181, 12.6560),
    "Billund": (55.7403, 9.1518),
    "Aarhus": (56.3000, 10.6190),
    "Aalborg": (57.0928, 9.8492),
    "Karup": (56.2975, 9.1246),
    "Esbjerg": (55.5259, 8.5534),
    "Bornholm": (55.0633, 14.7596),
    "Sønderborg": (54.9644, 9.7917),
    "Roskilde": (55.5856, 12.1314),
    "Thisted": (57.0688, 8.7052)
}

# -------------------------
# 3 Rounded curve function
# -------------------------

def rounded_route(start, end, bend=0.25, n_points=40):

    lat1, lon1 = start
    lat2, lon2 = end

    dx = lon2 - lon1
    dy = lat2 - lat1

    mid_lat = (lat1 + lat2) / 2
    mid_lon = (lon1 + lon2) / 2

    offset_lat = -dy * bend
    offset_lon = dx * bend

    ctrl_lat = mid_lat + offset_lat
    ctrl_lon = mid_lon + offset_lon

    points = []

    for t in np.linspace(0,1,n_points):

        lat = (1-t)**2 * lat1 + 2*(1-t)*t*ctrl_lat + t**2 * lat2
        lon = (1-t)**2 * lon1 + 2*(1-t)*t*ctrl_lon + t**2 * lon2

        points.append((lat,lon))

    return points


# -------------------------
# 4 Log scaling
# -------------------------

flows = df.values.flatten()
flows = flows[flows > 0]

min_log = np.log1p(flows.min())
max_log = np.log1p(flows.max())

def normalize(x):
    return (np.log1p(x) - min_log) / (max_log - min_log)


# -------------------------
# 5 Create map
# -------------------------

m = folium.Map(location=[56.2, 10.0], zoom_start=7)

# -------------------------
# 6 Draw airports
# -------------------------

for airport, coord in airport_coords.items():

    folium.CircleMarker(
        location=coord,
        radius=6,
        color="black",
        fill=True,
        fill_color="white",
        fill_opacity=1,
        popup=airport
    ).add_to(m)

# -------------------------
# 7 Red color scale
# -------------------------

colormap = cm.linear.Reds_09.scale(0,1)

# -------------------------
# 8 Draw OD routes
# -------------------------

for origin in df.index:

    for destination in df.columns:

        value = df.loc[origin, destination]

        if value <= 0:
            continue

        if origin == destination:
            continue

        start = airport_coords.get(origin)
        end = airport_coords.get(destination)

        if start is None or end is None:
            continue

        norm_val = normalize(value)

        curve = rounded_route(start, end)

        folium.PolyLine(
            curve,
            color=colormap(norm_val),
            weight=2 + 3*norm_val,
            opacity=0.9,
            dash_array="5,8",
            popup=f"{origin} → {destination}: {int(value):,}"
        ).add_to(m)

# -------------------------
# 9 Add legend
# -------------------------

colormap.caption = "Passenger flow (log scale)"
colormap.add_to(m)

# -------------------------
# 10 Save map
# -------------------------

m.save(BASE_DIR / "denmark_air_od_network_logscale.html")

print("Map saved as denmark_air_od_network_logscale.html")

# %% [markdown]
# 2020-2025

# %%
import pandas as pd
import folium
import numpy as np
import branca.colormap as cm

# -------------------------
# 1 Load Excel
# -------------------------

raw = pd.read_excel(BASE_DIR / "DKStatFly_years.xlsx", header=None)

# Destination names
destinations = raw.iloc[2,2:].str.replace("Til ", "", regex=False).tolist()

current_origin = None
records = []

# -------------------------
# 2 Parse Excel structure
# -------------------------

for i in range(3, len(raw)):

    origin_cell = raw.iloc[i,0]
    year = raw.iloc[i,1]

    # detect new origin
    if isinstance(origin_cell, str) and origin_cell.startswith("Fra"):
        current_origin = origin_cell.replace("Fra ", "")
        continue

    # detect year rows
    if pd.notna(year):

        for j, dest in enumerate(destinations):

            value = raw.iloc[i, j+2]

            records.append({
                "origin": current_origin,
                "destination": dest,
                "value": value
            })

# Convert to dataframe
df_long = pd.DataFrame(records)

# -------------------------
# 3 Build OD matrix
# -------------------------

od_matrix = (
    df_long
    .groupby(["origin","destination"])["value"]
    .sum()
    .unstack()
)

# remove "øvrige lufthavne"
if "øvrige lufthavne" in od_matrix.columns:
    od_matrix = od_matrix.drop(columns=["øvrige lufthavne"])

# clean string placeholders
od_matrix = od_matrix.replace("..........", np.nan)

# convert to numeric
od_matrix = od_matrix.apply(pd.to_numeric, errors="coerce")

# replace NaN with 0
od_matrix = od_matrix.fillna(0)

print("OD matrix used for map:")
print(od_matrix)

# -------------------------
# 4 Airport coordinates
# -------------------------

airport_coords = {
    "København": (55.6181, 12.6560),
    "Billund": (55.7403, 9.1518),
    "Aarhus": (56.3000, 10.6190),
    "Aalborg": (57.0928, 9.8492),
    "Karup": (56.2975, 9.1246),
    "Esbjerg": (55.5259, 8.5534),
    "Bornholm": (55.0633, 14.7596),
    "Sønderborg": (54.9644, 9.7917),
    "Roskilde": (55.5856, 12.1314),
    "Thisted": (57.0688, 8.7052)
}

# -------------------------
# 5 Rounded curve function
# -------------------------

def rounded_route(start, end, bend=0.25, n_points=40):

    lat1, lon1 = start
    lat2, lon2 = end

    dx = lon2 - lon1
    dy = lat2 - lat1

    mid_lat = (lat1 + lat2) / 2
    mid_lon = (lon1 + lon2) / 2

    offset_lat = -dy * bend
    offset_lon = dx * bend

    ctrl_lat = mid_lat + offset_lat
    ctrl_lon = mid_lon + offset_lon

    points = []

    for t in np.linspace(0,1,n_points):

        lat = (1-t)**2 * lat1 + 2*(1-t)*t*ctrl_lat + t**2 * lat2
        lon = (1-t)**2 * lon1 + 2*(1-t)*t*ctrl_lon + t**2 * lon2

        points.append((lat,lon))

    return points


# -------------------------
# 6 Log scaling
# -------------------------

flows = od_matrix.to_numpy().astype(float).flatten()
flows = flows[flows > 0]

min_log = np.log1p(flows.min())
max_log = np.log1p(flows.max())

def normalize(x):
    return (np.log1p(x) - min_log) / (max_log - min_log)

# -------------------------
# 7 Create map
# -------------------------

m = folium.Map(location=[56.2,10.0], zoom_start=7)

# Draw airports
for airport, coord in airport_coords.items():

    folium.CircleMarker(
        location=coord,
        radius=6,
        color="black",
        fill=True,
        fill_color="white",
        fill_opacity=1,
        popup=airport
    ).add_to(m)

# color scale
colormap = cm.linear.Reds_09.scale(0,1)

# -------------------------
# 8 Draw routes
# -------------------------

for origin in od_matrix.index:

    for destination in od_matrix.columns:

        value = od_matrix.loc[origin,destination]

        if value <= 0 or origin == destination:
            continue

        start = airport_coords.get(origin)
        end = airport_coords.get(destination)

        if start is None or end is None:
            continue

        norm_val = normalize(value)

        curve = rounded_route(start,end)

        folium.PolyLine(
            curve,
            color=colormap(norm_val),
            weight=2 + 3*norm_val,
            opacity=0.9,
            dash_array="5,8",
            popup=f"{origin} → {destination}: {int(value):,}"
        ).add_to(m)

# -------------------------
# 9 Add legend
# -------------------------

colormap.caption = "Passenger flow (sum over all years, log scale)"
colormap.add_to(m)

# -------------------------
# 10 Save map
# -------------------------

m.save(BASE_DIR / "denmark_air_od_network_all_years.html")

print("Map saved as denmark_air_od_network_all_years.html")

# %%
import pandas as pd
import folium
import numpy as np
import branca.colormap as cm

# -------------------------
# 1 Load Excel
# -------------------------

raw = pd.read_excel(BASE_DIR / "DKStatFly_years.xlsx", header=None)

# Destination names
destinations = raw.iloc[2,2:].str.replace("Til ", "", regex=False).tolist()

current_origin = None
records = []

# -------------------------
# 2 Parse Excel structure
# -------------------------

for i in range(3, len(raw)):

    origin_cell = raw.iloc[i,0]
    year = raw.iloc[i,1]

    # detect new origin
    if isinstance(origin_cell, str) and origin_cell.startswith("Fra"):

        current_origin = origin_cell.replace("Fra ", "")

        # IMPORTANT: read the values on the SAME ROW (this captures 2020)
        for j, dest in enumerate(destinations):

            value = raw.iloc[i, j+2]

            records.append({
                "origin": current_origin,
                "destination": dest,
                "value": value
            })

        continue

    # detect year rows
    if pd.notna(year):

        for j, dest in enumerate(destinations):

            value = raw.iloc[i, j+2]

            records.append({
                "origin": current_origin,
                "destination": dest,
                "value": value
            })

# Convert to dataframe
df_long = pd.DataFrame(records)

# -------------------------
# 3 Build OD matrix
# -------------------------

od_matrix = (
    df_long
    .groupby(["origin","destination"])["value"]
    .sum()
    .unstack()
)

# remove "øvrige lufthavne"
if "øvrige lufthavne" in od_matrix.columns:
    od_matrix = od_matrix.drop(columns=["øvrige lufthavne"])

# clean string placeholders
od_matrix = od_matrix.replace("..........", np.nan)
od_matrix = od_matrix.replace("..", np.nan)

# convert to numeric
od_matrix = od_matrix.apply(pd.to_numeric, errors="coerce")

# replace NaN with 0
od_matrix = od_matrix.fillna(0)

# -------------------------
# NEW: convert thousands → passengers
# -------------------------

od_matrix = od_matrix * 1000

print("OD matrix used for map:")
print(od_matrix)

# -------------------------
# 4 Airport coordinates
# -------------------------

airport_coords = {
    "København": (55.6181, 12.6560),
    "Billund": (55.7403, 9.1518),
    "Aarhus": (56.3000, 10.6190),
    "Aalborg": (57.0928, 9.8492),
    "Karup": (56.2975, 9.1246),
    "Esbjerg": (55.5259, 8.5534),
    "Bornholm": (55.0633, 14.7596),
    "Sønderborg": (54.9644, 9.7917),
    "Roskilde": (55.5856, 12.1314),
    "Thisted": (57.0688, 8.7052)
}

# -------------------------
# 5 Rounded curve function
# -------------------------

def rounded_route(start, end, bend=0.25, n_points=40):

    lat1, lon1 = start
    lat2, lon2 = end

    dx = lon2 - lon1
    dy = lat2 - lat1

    mid_lat = (lat1 + lat2) / 2
    mid_lon = (lon1 + lon2) / 2

    offset_lat = -dy * bend
    offset_lon = dx * bend

    ctrl_lat = mid_lat + offset_lat
    ctrl_lon = mid_lon + offset_lon

    points = []

    for t in np.linspace(0,1,n_points):

        lat = (1-t)**2 * lat1 + 2*(1-t)*t*ctrl_lat + t**2 * lat2
        lon = (1-t)**2 * lon1 + 2*(1-t)*t*ctrl_lon + t**2 * lon2

        points.append((lat,lon))

    return points


# -------------------------
# 6 Log scaling
# -------------------------

flows = od_matrix.to_numpy().astype(float).flatten()
flows = flows[flows > 0]

min_log = np.log1p(flows.min())
max_log = np.log1p(flows.max())

def normalize(x):
    return (np.log1p(x) - min_log) / (max_log - min_log)

# -------------------------
# 7 Create map
# -------------------------

m = folium.Map(location=[56.2,10.0], zoom_start=7)

# Draw airports
for airport, coord in airport_coords.items():

    folium.CircleMarker(
        location=coord,
        radius=6,
        color="black",
        fill=True,
        fill_color="white",
        fill_opacity=1,
        popup=airport
    ).add_to(m)

# color scale
colormap = cm.linear.Reds_09.scale(0,1)

# -------------------------
# 8 Draw routes
# -------------------------

for origin in od_matrix.index:

    for destination in od_matrix.columns:

        value = od_matrix.loc[origin,destination]

        if value <= 0 or origin == destination:
            continue

        start = airport_coords.get(origin)
        end = airport_coords.get(destination)

        if start is None or end is None:
            continue

        norm_val = normalize(value)

        curve = rounded_route(start,end)

        folium.PolyLine(
            curve,
            color=colormap(norm_val),
            weight=2 + 3*norm_val,
            opacity=0.9,
            dash_array="5,8",
            popup=f"{origin} → {destination}: {int(value):,}"
        ).add_to(m)

# -------------------------
# 9 Add legend
# -------------------------

colormap.caption = "Passenger flow (sum over all years, log scale)"
colormap.add_to(m)

# -------------------------
# 10 Save map
# -------------------------

m.save(BASE_DIR / "denmark_air_od_network_all_years.html")

print("Map saved as denmark_air_od_network_all_years.html")


