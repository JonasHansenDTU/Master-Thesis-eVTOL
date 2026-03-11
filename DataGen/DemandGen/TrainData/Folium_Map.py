import json
from pathlib import Path
import re
import unicodedata

import folium
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
OD_CSV = BASE_DIR / 'OD_filtered_pruned.csv'
OD_XLSX = BASE_DIR / 'OD_filtered_pruned.xlsx'
EXPORT_JSON = BASE_DIR / 'export.json'
OUT_HTML = BASE_DIR / 'OD_filtered_pruned_map.html'

MIN_OD_VALUE = 1
TOP_LINES_PER_LAYER = 400
THRESHOLD_LAYERS = [500, 2000, 10000]
KEEP_STRONGER_DIRECTION_ONLY = True

# Optional origin-focus mode (set station name, or None for full network)
FOCUS_ORIGIN = None
FOCUS_INCLUDE_INBOUND = True


def normalize_station_name(name: str) -> str:
	text = str(name).strip().lower()
	text = unicodedata.normalize('NFKD', text)
	text = ''.join(ch for ch in text if not unicodedata.combining(ch))
	text = re.sub(r'\bst\.?\b', 'station', text)
	text = re.sub(r'\s+', ' ', text)
	return text


def load_od_matrix() -> pd.DataFrame:
	if OD_CSV.exists():
		od_df = pd.read_csv(OD_CSV, index_col=0)
	elif OD_XLSX.exists():
		od_df = pd.read_excel(OD_XLSX, index_col=0)
	else:
		raise FileNotFoundError('Could not find OD_filtered_pruned.csv or OD_filtered_pruned.xlsx')

	od_df = od_df.apply(pd.to_numeric, errors='coerce').fillna(0)
	od_df = od_df.loc[od_df.index, od_df.index.intersection(od_df.columns)]
	od_df = od_df.reindex(index=od_df.columns, columns=od_df.columns, fill_value=0)
	return od_df


def load_station_coordinates() -> dict[str, tuple[float, float]]:
	with open(EXPORT_JSON, 'r', encoding='utf-8') as handle:
		osm_data = json.load(handle)

	coord_lookup: dict[str, tuple[float, float]] = {}
	for element in osm_data.get('elements', []):
		tags = element.get('tags', {})
		name = tags.get('name')
		lat = element.get('lat')
		lon = element.get('lon')
		if name is None or lat is None or lon is None:
			continue
		norm_name = normalize_station_name(name)
		if norm_name and norm_name not in coord_lookup:
			coord_lookup[norm_name] = (float(lat), float(lon))

	return coord_lookup


def build_edges(od_matrix: pd.DataFrame, stations_map: list[str]) -> list[tuple[str, str, float]]:
	edges_directed: list[tuple[str, str, float]] = []

	focus_origin = FOCUS_ORIGIN if FOCUS_ORIGIN in stations_map else None

	for origin in stations_map:
		for destination in stations_map:
			if origin == destination:
				continue

			value = float(od_matrix.at[origin, destination])
			if value < MIN_OD_VALUE:
				continue

			if focus_origin is not None:
				if origin == focus_origin:
					edges_directed.append((origin, destination, value))
				elif FOCUS_INCLUDE_INBOUND and destination == focus_origin:
					edges_directed.append((origin, destination, value))
			else:
				edges_directed.append((origin, destination, value))

	if KEEP_STRONGER_DIRECTION_ONLY:
		pair_best: dict[tuple[str, str], tuple[str, str, float]] = {}
		for origin, destination, value in edges_directed:
			key = tuple(sorted((origin, destination)))
			current = pair_best.get(key)
			if current is None or value > current[2]:
				pair_best[key] = (origin, destination, value)
		edges = list(pair_best.values())
	else:
		edges = edges_directed

	edges.sort(key=lambda item: item[2], reverse=True)
	return edges


def add_threshold_layers(
	od_map: folium.Map,
	edges: list[tuple[str, str, float]],
	coords_by_station: dict[str, tuple[float, float]],
):
	for threshold in THRESHOLD_LAYERS:
		layer_edges = [edge for edge in edges if edge[2] >= threshold][:TOP_LINES_PER_LAYER]
		group = folium.FeatureGroup(name=f'OD ≥ {threshold:,} ({len(layer_edges)} lines)', show=(threshold == THRESHOLD_LAYERS[0]))

		if layer_edges:
			values = np.array([edge[2] for edge in layer_edges], dtype=float)
			log_values = np.log1p(values)
			vmin = float(log_values.min())
			vmax = float(log_values.max())

			for (origin, destination, value), log_value in zip(layer_edges, log_values):
				lat1, lon1 = coords_by_station[origin]
				lat2, lon2 = coords_by_station[destination]

				if vmax > vmin:
					norm = (log_value - vmin) / (vmax - vmin)
				else:
					norm = 0.5

				weight = 1.0 + 6.0 * norm
				opacity = 0.15 + 0.55 * norm

				folium.PolyLine(
					locations=[[lat1, lon1], [lat2, lon2]],
					weight=weight,
					opacity=opacity,
					color='#1f78b4',
					tooltip=f'{origin} → {destination}: {value:,.0f}',
				).add_to(group)

		group.add_to(od_map)


def main():
	od_matrix = load_od_matrix()
	coord_lookup = load_station_coordinates()

	stations = list(od_matrix.index)
	coords_by_station: dict[str, tuple[float, float]] = {}
	for station in stations:
		norm_station = normalize_station_name(station)
		if norm_station in coord_lookup:
			coords_by_station[station] = coord_lookup[norm_station]

	stations_map = [station for station in stations if station in coords_by_station]
	if not stations_map:
		raise ValueError('No stations from OD_filtered_pruned matched coordinates in export.json')

	center_lat = float(np.mean([coords_by_station[s][0] for s in stations_map]))
	center_lon = float(np.mean([coords_by_station[s][1] for s in stations_map]))

	od_map = folium.Map(location=[center_lat, center_lon], zoom_start=7, tiles='CartoDB positron')

	station_group = folium.FeatureGroup(name='Stations', show=True)
	for station in stations_map:
		lat, lon = coords_by_station[station]
		folium.CircleMarker(
			location=[lat, lon],
			radius=2,
			color='black',
			fill=True,
			fill_opacity=0.9,
			popup=station,
		).add_to(station_group)
	station_group.add_to(od_map)

	edges = build_edges(od_matrix, stations_map)
	add_threshold_layers(od_map, edges, coords_by_station)
	folium.LayerControl(collapsed=False).add_to(od_map)

	od_map.save(str(OUT_HTML))
	print(f'Map saved to: {OUT_HTML.name}')
	print(f'Stations plotted: {len(stations_map)} / {len(stations)}')
	print(f'Candidate OD lines after filtering: {len(edges)}')
	print(f'Layers: {THRESHOLD_LAYERS} with top {TOP_LINES_PER_LAYER} lines per layer')
	if FOCUS_ORIGIN is not None:
		print(f'Focus origin mode: {FOCUS_ORIGIN} (inbound included: {FOCUS_INCLUDE_INBOUND})')


if __name__ == '__main__':
	main()
