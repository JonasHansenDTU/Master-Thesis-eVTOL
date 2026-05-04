# %%
import pandas as pd
import json

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
# 2 GeoJSON navne
# -------------------------

features = [
    f for f in geo["features"]
    if f["geometry"]["type"] in ["Polygon", "MultiPolygon"]
]

geo_names = set(
    f["properties"].get("name")
    for f in features
    if f["properties"].get("name") is not None
)

# -------------------------
# 3 Lav Geo-navne fra data
# -------------------------

koder["navn_geo"] = koder["navn"] + " Kommune"

# Special cases
koder.loc[koder["navn"] == "København", "navn_geo"] = "Københavns Kommune"
koder.loc[koder["navn"] == "Bornholm", "navn_geo"] = "Bornholms Regionskommune"

# Fjern Christiansø (findes ikke i geojson)
koder = koder[koder["navn"] != "Christiansø"]

# -------------------------
# 4 Mapping (kode → geo navn)
# -------------------------

kode_to_geo = dict(zip(koder["kode"], koder["navn_geo"]))

df["FraNavn"] = df["FraKommune"].map(kode_to_geo)
df["TilNavn"] = df["TilKommune"].map(kode_to_geo)

# -------------------------
# 5 Tjek for manglende mapping
# -------------------------

missing_fra = df[df["FraNavn"].isna()]["FraKommune"].unique()
missing_til = df[df["TilNavn"].isna()]["TilKommune"].unique()

print("\n--- MANGLENDE KODER ---")
print("FraKommune:", missing_fra)
print("TilKommune:", missing_til)

# -------------------------
# 6 Sammenlign navne
# -------------------------

data_names = set(koder["navn_geo"])

missing_in_geo = data_names - geo_names
missing_in_data = geo_names - data_names

print("\n--- MATCH CHECK ---")
print("Mangler i GeoJSON:")
print(missing_in_geo)

print("\nGeoJSON uden match i data:")
print(missing_in_data)

# -------------------------
# 7 Final check
# -------------------------

valid_rows = df[
    df["FraNavn"].isin(geo_names) &
    df["TilNavn"].isin(geo_names)
]

print("\n--- FINAL STATUS ---")
print(f"Total rækker: {len(df)}")
print(f"Gyldige rækker: {len(valid_rows)}")

if len(valid_rows) == len(df):
    print("ALT MATCHER PERFEKT")
else:
    print("Der er stadig mismatch")

print("\nDone ")