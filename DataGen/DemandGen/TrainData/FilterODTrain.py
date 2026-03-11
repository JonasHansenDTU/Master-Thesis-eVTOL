import pandas as pd
import matplotlib.pyplot as plt
import json
import numpy as np
import re
import unicodedata
from pathlib import Path


def save_matrix(df_obj: pd.DataFrame, excel_name: str, csv_name: str):
	df_obj.to_csv(csv_name, encoding='utf-8-sig')
	try:
		df_obj.to_excel(excel_name)
	except PermissionError:
		excel_path = Path(excel_name)
		fallback = excel_path.with_stem(f'{excel_path.stem}_new')
		df_obj.to_excel(fallback)
		print(f'File locked: {excel_name}. Saved Excel to {fallback.name} instead.')


def normalize_station_name(name: str) -> str:
	if pd.isna(name):
		return ''
	text = str(name).strip().lower()
	text = unicodedata.normalize('NFKD', text)
	text = ''.join(ch for ch in text if not unicodedata.combining(ch))
	text = re.sub(r'\bst\.?\b', 'station', text)
	text = re.sub(r'\s+', ' ', text)
	return text


def haversine_km(lat1, lon1, lat2, lon2):
	R = 6371.0
	phi1 = np.radians(lat1)
	phi2 = np.radians(lat2)
	dphi = np.radians(lat2 - lat1)
	dlambda = np.radians(lon2 - lon1)
	a = np.sin(dphi / 2.0) ** 2 + np.cos(phi1) * np.cos(phi2) * np.sin(dlambda / 2.0) ** 2
	c = 2 * np.arctan2(np.sqrt(a), np.sqrt(1 - a))
	return R * c


# Load OD data
df = pd.read_excel('DataGen/DemandGen/TrainData/OD_eksport_samlet.xlsx')
df.columns = df.columns.str.strip()

from_col = next(c for c in df.columns if c.lower().startswith('fra station'))
to_col = next(c for c in df.columns if c.lower().startswith('til station'))
trips_col = next(c for c in df.columns if c.lower().startswith('antal rejser'))

od_matrix = pd.pivot_table(
	df,
	index=from_col,
	columns=to_col,
	values=trips_col,
	aggfunc='sum',
	fill_value=0,
)
od_matrix = od_matrix.sort_index().sort_index(axis=1)

print(od_matrix)
save_matrix(od_matrix, 'DataGen/DemandGen/TrainData/OD_matrix.xlsx', 'DataGen/DemandGen/TrainData/OD_matrix.csv')


# # Optional: Simple heat map of OD matrix
# plt.figure(figsize=(10, 8))
# plt.imshow(od_matrix.values, aspect='auto', cmap='viridis', interpolation='nearest')
# plt.colorbar(label='Antal Rejser')
# plt.title('OD Heatmap (Station to Station)')
# plt.xlabel('Til Station')
# plt.ylabel('Fra Station')
# plt.xticks([])
# plt.yticks([])
# plt.tight_layout()
# plt.savefig('OD_heatmap.png', dpi=200)
# plt.show()


# Build station distance matrix from export.json
with open('DataGen/DemandGen/TrainData/export.json', 'r', encoding='utf-8') as f:
	osm_data = json.load(f)

coord_lookup = {}
for element in osm_data.get('elements', []):
	tags = element.get('tags', {})
	name = tags.get('name')
	lat = element.get('lat')
	lon = element.get('lon')
	if name is None or lat is None or lon is None:
		continue
	norm_name = normalize_station_name(name)
	if norm_name and norm_name not in coord_lookup:
		coord_lookup[norm_name] = (name, lat, lon)

stations = sorted(set(od_matrix.index).union(set(od_matrix.columns)))

coords_by_station = {}
for station in stations:
	norm_station = normalize_station_name(station)
	if norm_station in coord_lookup:
		_, lat, lon = coord_lookup[norm_station]
		coords_by_station[station] = (lat, lon)

distance_matrix = pd.DataFrame(np.nan, index=stations, columns=stations, dtype=float)
for i, station_i in enumerate(stations):
	if station_i not in coords_by_station:
		continue
	lat1, lon1 = coords_by_station[station_i]
	for j, station_j in enumerate(stations):
		if station_j not in coords_by_station:
			continue
		if i == j:
			distance_matrix.iat[i, j] = 0.0
			continue
		lat2, lon2 = coords_by_station[station_j]
		distance_matrix.iat[i, j] = haversine_km(lat1, lon1, lat2, lon2)

save_matrix(distance_matrix, 'DataGen/DemandGen/TrainData/Station_distance_matrix_km.xlsx', 'DataGen/DemandGen/TrainData/Station_distance_matrix_km.csv')

matched_count = len(coords_by_station)
print(f'Matched stations with coordinates: {matched_count}/{len(stations)}')

missing_stations = [s for s in stations if s not in coords_by_station]
if missing_stations:
	print('Stations missing coordinates (first 20):')
	print(missing_stations[:20])


# Build NEW filtered OD matrix only (original od_matrix is unchanged)
DISTANCE_THRESHOLD_KM = 50.0  # Change X here

stations_with_coords = [s for s in stations if s in coords_by_station]
common_stations = [s for s in stations_with_coords if s in od_matrix.index and s in od_matrix.columns]

od_matrix_subset = od_matrix.reindex(index=common_stations, columns=common_stations, fill_value=0).copy()
distance_submatrix = distance_matrix.reindex(index=common_stations, columns=common_stations)
distance_submatrix = distance_submatrix.apply(pd.to_numeric, errors='coerce')

close_mask = distance_submatrix.lt(DISTANCE_THRESHOLD_KM) & distance_submatrix.notna()

OD_filtered = od_matrix_subset.mask(close_mask, -1)

# Further filtering: keep only stations i with at least one j where OD_filtered[i, j] > threshold
OD_VALUE_THRESHOLD = OD_filtered.values.max() * 0.2  # 20% of maximum value

keep_mask = OD_filtered.gt(OD_VALUE_THRESHOLD).any(axis=1)
stations_to_keep = OD_filtered.index[keep_mask]

# Keep matrix square (same station set on rows/columns)
OD_filtered_pruned = OD_filtered.loc[stations_to_keep, stations_to_keep].copy()

save_matrix(OD_filtered, 'DataGen/DemandGen/TrainData/OD_filtered.xlsx', 'DataGen/DemandGen/TrainData/OD_filtered.csv')
save_matrix(OD_filtered_pruned, 'DataGen/DemandGen/TrainData/OD_filtered_pruned.xlsx', 'DataGen/DemandGen/TrainData/OD_filtered_pruned.csv')

print(f'New OD matrix shape: {OD_filtered.shape}')
print(f'Pairs set to -1 (distance < {DISTANCE_THRESHOLD_KM} km): {int(close_mask.sum().sum())}')
print(f'Pruned OD matrix shape (OD > {OD_VALUE_THRESHOLD} exists in row): {OD_filtered_pruned.shape}')
print('Original od_matrix is unchanged. Close-station edits are only in OD_filtered.')