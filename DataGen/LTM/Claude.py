"""
Danish Municipality OD Flow Map
================================
Visualizes origin-destination flows between Danish municipalities
for four transport modes: bil, kollektiv, varebil, and fly.

Requirements:
    pip install pandas geopandas folium openpyxl shapely

Usage:
    1. Set EXCEL_PATH and GEOJSON_PATH below.
    2. Run: python od_flow_map.py
    3. Open od_flow_map.html in a browser.
"""

import json
import numpy as np
import pandas as pd
import geopandas as gpd
import folium
from folium import plugins
from pathlib import Path


# ─────────────────────────────────────────────
# CONFIG  —  edit these paths
# ─────────────────────────────────────────────
EXCEL_PATH = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"
GEOJSON_PATH = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/kommune.geojson"



OUTPUT_HTML  = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/od_flow_map.html"

# Column names in the Excel file
COL_FROM = "FraKommune"
COL_TO   = "TilKommune"

# GeoJSON property that holds the municipality name
GEOJSON_NAME_PROP = "name"

# Map each municipality code (int) → GeoJSON name string.
# Extend / replace this dict to match your actual data.
# Example:  185 → "Frederiksberg Kommune"
# If you already have a full lookup table as a CSV/sheet, load it instead.
KOMMUNE_CODE_TO_NAME: dict[int, str] = {
    # ── Greater Copenhagen ──────────────────────────────────────────────
    101: "Københavns Kommune",
    147: "Frederiksberg Kommune",
    151: "Ballerup Kommune",
    153: "Brøndby Kommune",
    155: "Dragør Kommune",
    157: "Gentofte Kommune",
    159: "Gladsaxe Kommune",
    161: "Glostrup Kommune",
    163: "Herlev Kommune",
    165: "Albertslund Kommune",
    167: "Hvidovre Kommune",
    169: "Høje-Taastrup Kommune",
    173: "Lyngby-Taarbæk Kommune",
    175: "Rødovre Kommune",
    183: "Ishøj Kommune",
    185: "Tårnby Kommune",          # CPH airport municipality
    187: "Vallensbæk Kommune",
    190: "Furesø Kommune",
    201: "Allerød Kommune",
    210: "Fredensborg Kommune",
    217: "Helsingør Kommune",
    219: "Hillerød Kommune",
    223: "Hørsholm Kommune",
    230: "Rudersdal Kommune",
    240: "Egedal Kommune",
    250: "Frederikssund Kommune",
    253: "Greve Kommune",
    259: "Køge Kommune",
    260: "Halsnæs Kommune",
    265: "Roskilde Kommune",
    269: "Solrød Kommune",
    270: "Gribskov Kommune",
    # ── Bornholm ────────────────────────────────────────────────────────
    400: "Bornholms Regionskommune",
    # ── South Denmark ───────────────────────────────────────────────────
    410: "Middelfart Kommune",
    420: "Assens Kommune",
    430: "Faaborg-Midtfyn Kommune",
    440: "Kerteminde Kommune",
    450: "Nyborg Kommune",
    461: "Odense Kommune",
    479: "Svendborg Kommune",
    480: "Nordfyns Kommune",
    482: "Langeland Kommune",
    492: "Ærø Kommune",
    510: "Haderslev Kommune",
    530: "Billund Kommune",         # BLL airport municipality
    540: "Sønderborg Kommune",      # SGD airport municipality
    550: "Tønder Kommune",
    561: "Esbjerg Kommune",
    563: "Fanø Kommune",
    573: "Varde Kommune",
    575: "Vejen Kommune",
    580: "Aabenraa Kommune",
    607: "Fredericia Kommune",
    615: "Horsens Kommune",
    621: "Kolding Kommune",
    630: "Vejle Kommune",
    # ── Central Jutland ─────────────────────────────────────────────────
    657: "Herning Kommune",
    661: "Holstebro Kommune",
    665: "Lemvig Kommune",
    671: "Struer Kommune",
    706: "Syddjurs Kommune",        # Near AAR airport
    707: "Norddjurs Kommune",
    710: "Favrskov Kommune",
    727: "Odder Kommune",
    730: "Randers Kommune",
    740: "Silkeborg Kommune",
    741: "Samsø Kommune",
    746: "Skanderborg Kommune",
    751: "Aarhus Kommune",
    756: "Ikast-Brande Kommune",
    760: "Ringkøbing-Skjern Kommune",
    766: "Hedensted Kommune",
    773: "Morsø Kommune",
    779: "Skive Kommune",
    787: "Thisted Kommune",
    791: "Viborg Kommune",          # BLL catchment
    # ── North Jutland ───────────────────────────────────────────────────
    810: "Brønderslev Kommune",
    813: "Frederikshavn Kommune",
    820: "Vesthimmerlands Kommune",
    825: "Læsø Kommune",
    840: "Rebild Kommune",
    846: "Mariagerfjord Kommune",
    849: "Jammerbugt Kommune",
    851: "Aalborg Kommune",         # AAL airport municipality
    860: "Hjørring Kommune",
    # ── Zealand ─────────────────────────────────────────────────────────
    306: "Odsherred Kommune",
    316: "Holbæk Kommune",
    320: "Faxe Kommune",
    326: "Kalundborg Kommune",
    329: "Ringsted Kommune",
    330: "Slagelse Kommune",
    336: "Stevns Kommune",
    340: "Sorø Kommune",
    350: "Lejre Kommune",
    360: "Lolland Kommune",
    370: "Næstved Kommune",
    376: "Guldborgsund Kommune",
    390: "Vordingborg Kommune",
}

# Transport modes to visualise
MODES = {
    "bil":       {"color": "#FF6B35", "label": "Bil",        "top_n": 60},
    "kollektiv": {"color": "#4ECDC4", "label": "Kollektiv",  "top_n": 60},
    "varebil":   {"color": "#FFE66D", "label": "Varebil",    "top_n": 40},
    "fly":       {"color": "#00BFFF", "label": "Fly (alle)", "top_n": 9999},
}


# ─────────────────────────────────────────────
# HELPERS
# ─────────────────────────────────────────────

def bezier_arc(
    p1: tuple[float, float],
    p2: tuple[float, float],
    n: int = 60,
    curvature: float = 0.25,
) -> list[tuple[float, float]]:
    """Quadratic Bézier arc between two (lat, lon) points."""
    lat1, lon1 = p1
    lat2, lon2 = p2
    mid_lat = (lat1 + lat2) / 2
    mid_lon = (lon1 + lon2) / 2
    dlat = lat2 - lat1
    dlon = lon2 - lon1
    ctrl_lat = mid_lat - curvature * dlon
    ctrl_lon = mid_lon + curvature * dlat
    t = np.linspace(0, 1, n)
    lats = (1 - t) ** 2 * lat1 + 2 * (1 - t) * t * ctrl_lat + t**2 * lat2
    lons = (1 - t) ** 2 * lon1 + 2 * (1 - t) * t * ctrl_lon + t**2 * lon2
    return list(zip(lats.tolist(), lons.tolist()))


def scale_width(value: float, max_value: float, min_px: float = 1, max_px: float = 8) -> float:
    """Scale a flow value to a line-width in pixels."""
    if max_value == 0:
        return min_px
    return min_px + (value / max_value) * (max_px - min_px)


def run_diagnostics(df: pd.DataFrame) -> None:
    """Print a quick diagnostic summary for every mode column."""
    print("\n" + "=" * 60)
    print("DIAGNOSTICS")
    print("=" * 60)
    for mode in MODES:
        col = f"AntalErhvervsture_{mode}"
        if col not in df.columns:
            print(f"  [{mode}]  ⚠  column '{col}' not found in data")
            continue
        nz = (df[col] > 0).sum()
        total = df[col].sum()
        print(f"\n  [{mode}]  non-zero rows: {nz}/{len(df)}   sum: {total:,.1f}")
        if nz > 0:
            top = (
                df[df[col] > 0]
                .nlargest(5, col)[[COL_FROM, COL_TO, col]]
                .to_string(index=False)
            )
            print(top)
    print("=" * 60 + "\n")


# ─────────────────────────────────────────────
# DATA LOADING
# ─────────────────────────────────────────────

def load_data(excel_path: str, geojson_path: str):
    """Load Excel OD data and GeoJSON boundaries; compute centroids."""

    # ── Excel ──────────────────────────────────────────────────────────
    print(f"Loading Excel:   {excel_path}")
    df = pd.read_excel(excel_path)
    print(f"  Rows: {len(df):,}   Columns: {list(df.columns)}")

    # Ensure code columns are integers
    df[COL_FROM] = df[COL_FROM].astype(int)
    df[COL_TO]   = df[COL_TO].astype(int)

    # ── GeoJSON ────────────────────────────────────────────────────────
    print(f"Loading GeoJSON: {geojson_path}")
    gdf = gpd.read_file(geojson_path)
    if gdf.crs and gdf.crs.to_epsg() != 4326:
        gdf = gdf.to_crs(epsg=4326)

    # Compute centroids (project to a metric CRS first for accuracy)
    gdf_proj = gdf.to_crs(epsg=25832)          # UTM zone 32N — good for Denmark
    gdf["centroid_lat"] = gdf_proj.geometry.centroid.to_crs(epsg=4326).y
    gdf["centroid_lon"] = gdf_proj.geometry.centroid.to_crs(epsg=4326).x

    # Build centroid lookup:  kommune_name → (lat, lon)
    centroids_by_name: dict[str, tuple[float, float]] = {
        row[GEOJSON_NAME_PROP]: (row["centroid_lat"], row["centroid_lon"])
        for _, row in gdf.iterrows()
        if row[GEOJSON_NAME_PROP]
    }

    # Build centroid lookup by code using KOMMUNE_CODE_TO_NAME
    centroids_by_code: dict[int, tuple[float, float]] = {}
    missing_codes: list[int] = []

    all_codes = set(df[COL_FROM].unique()) | set(df[COL_TO].unique())
    for code in all_codes:
        name = KOMMUNE_CODE_TO_NAME.get(code)
        if name and name in centroids_by_name:
            centroids_by_code[code] = centroids_by_name[name]
        else:
            missing_codes.append(code)

    if missing_codes:
        print(f"  ⚠  {len(missing_codes)} codes have no centroid mapping: {sorted(missing_codes)}")
        print("     Add them to KOMMUNE_CODE_TO_NAME at the top of this script.")

    return df, gdf, centroids_by_code


# ─────────────────────────────────────────────
# MAP BUILDING
# ─────────────────────────────────────────────

def build_map(
    df: pd.DataFrame,
    gdf: gpd.GeoDataFrame,
    centroids: dict[int, tuple[float, float]],
) -> folium.Map:

    # ── Base map ───────────────────────────────────────────────────────
    m = folium.Map(
        location=[56.0, 10.5],
        zoom_start=7,
        tiles="CartoDB voyager",
        prefer_canvas=True,
    )

    # ── Municipality boundaries ────────────────────────────────────────
    # Filter to polygon features only — this prevents Folium from auto-adding
    # a centroid pin marker for every feature when a tooltip is attached.
    geo_interface = gdf.__geo_interface__
    polygon_features = [
        f for f in geo_interface["features"]
        if f["geometry"]["type"] in ("Polygon", "MultiPolygon")
    ]
    polygon_collection = {"type": "FeatureCollection", "features": polygon_features}

    boundary_layer = folium.FeatureGroup(name="Kommunegrænser", show=True)
    folium.GeoJson(
        polygon_collection,
        style_function=lambda _: {
            "fillColor":   "#2a2a2a",
            "color":       "#666666",
            "weight":      1.0,
            "fillOpacity": 0.25,
        },
        highlight_function=lambda _: {
            "fillColor":   "#3a3a3a",
            "color":       "#bbbbbb",
            "weight":      2.0,
            "fillOpacity": 0.5,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=[GEOJSON_NAME_PROP],
            aliases=["Kommune:"],
            localize=True,
        ),
    ).add_to(boundary_layer)
    boundary_layer.add_to(m)

    # ── One FeatureGroup per transport mode ────────────────────────────
    for mode, cfg in MODES.items():
        col = f"AntalErhvervsture_{mode}"
        if col not in df.columns:
            print(f"  Skipping {mode}: column not found.")
            continue

        df_sub = df[df[col] > 0].copy()

        # Keep only rows where both endpoints have centroids
        df_sub = df_sub[
            df_sub[COL_FROM].isin(centroids) & df_sub[COL_TO].isin(centroids)
        ]

        if df_sub.empty:
            print(f"  [{mode}]  No plottable rows (all zeros or missing centroids).")
            continue

        # For dense modes, show only the top-N flows to avoid clutter
        if cfg["top_n"] < len(df_sub):
            df_sub = df_sub.nlargest(cfg["top_n"], col)

        max_flow = df_sub[col].max()
        fg = folium.FeatureGroup(name=f"{cfg['label']} ({len(df_sub)} flows)", show=False)

        for _, row in df_sub.iterrows():
            fra_coord = centroids[row[COL_FROM]]
            til_coord = centroids[row[COL_TO]]
            flow      = row[col]

            arc   = bezier_arc(fra_coord, til_coord, curvature=0.28)
            width = scale_width(flow, max_flow)

            fra_name = KOMMUNE_CODE_TO_NAME.get(row[COL_FROM], str(row[COL_FROM]))
            til_name = KOMMUNE_CODE_TO_NAME.get(row[COL_TO],   str(row[COL_TO]))

            folium.PolyLine(
                arc,
                weight=width,
                color=cfg["color"],
                opacity=0.65,
                tooltip=(
                    f"<b>{fra_name}</b> → <b>{til_name}</b><br>"
                    f"Erhvervsture ({mode}): <b>{flow:,.0f}</b>"
                ),
            ).add_to(fg)

        fg.add_to(m)

    # ── Show "bil" layer by default ────────────────────────────────────
    # Toggle it on in the LayerControl by rebuilding with show=True
    # (already set in MODES — bil is drawn first above)
    # Re-enable bil layer visibility:
    for child in m._children.values():
        if isinstance(child, folium.FeatureGroup) and child.layer_name.startswith("Bil"):
            child.show = True

    # ── Legend ─────────────────────────────────────────────────────────
    legend_html = """
    <div style="
        position: fixed; bottom: 40px; left: 40px; z-index: 1000;
        background: rgba(20,20,20,0.85); border-radius: 10px;
        padding: 14px 18px; color: white; font-family: sans-serif; font-size: 13px;
        box-shadow: 0 2px 8px rgba(0,0,0,0.5);
    ">
        <b style="font-size:14px;">Erhvervsture — transportmiddel</b><br><br>
        <span style="color:#FF6B35;">&#9644;</span> Bil<br>
        <span style="color:#4ECDC4;">&#9644;</span> Kollektiv<br>
        <span style="color:#FFE66D;">&#9644;</span> Varebil<br>
        <span style="color:#00BFFF;">&#9644;</span> Fly<br>
        <br>
        <span style="opacity:0.6; font-size:11px;">Linjetykkelse ~ antal ture</span>
    </div>
    """
    m.get_root().html.add_child(folium.Element(legend_html))

    # ── Layer control ──────────────────────────────────────────────────
    folium.LayerControl(collapsed=False).add_to(m)

    return m


# ─────────────────────────────────────────────
# MAIN
# ─────────────────────────────────────────────

def main() -> None:
    # Validate paths
    for path, label in [(EXCEL_PATH, "Excel"), (GEOJSON_PATH, "GeoJSON")]:
        if not Path(path).exists():
            raise FileNotFoundError(
                f"{label} file not found: '{path}'\n"
                "Update EXCEL_PATH / GEOJSON_PATH at the top of the script."
            )

    df, gdf, centroids = load_data(EXCEL_PATH, GEOJSON_PATH)

    run_diagnostics(df)

    print("Building map …")
    m = build_map(df, gdf, centroids)

    m.save(OUTPUT_HTML)
    print(f"\n✅  Map saved → {OUTPUT_HTML}")
    print("   Open it in any browser.")


if __name__ == "__main__":
    main()