# %%
import pandas as pd
import numpy as np
from pathlib import Path

# -------------------------
# Paths
# -------------------------

BASE_DIR = Path("/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM")
file_path = BASE_DIR / "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"

# -------------------------
# Load OD data
# -------------------------

df = pd.read_excel(file_path, sheet_name="LTM")

MODE = "fly"
col = f"AntalErhvervsture_{MODE}"

df = df[["FraKommune", "TilKommune", col]].copy()
df.rename(columns={col: "trips"}, inplace=True)

# remove zero flows
df = df[df["trips"] > 0]

# remove self trips
df = df[df["FraKommune"] != df["TilKommune"]]

# -------------------------
# Load vertiports
# -------------------------

df_verti = pd.read_excel(file_path, sheet_name="VertiportID")

kommune_to_vertiport = dict(zip(df_verti["Kommuneid"], df_verti["id"]))
all_vertiports = list(df_verti["id"])

# keep only routes with vertiports
df = df[
    df["FraKommune"].isin(kommune_to_vertiport) &
    df["TilKommune"].isin(kommune_to_vertiport)
]

print(f"Ruter efter vertiport filter: {len(df)}")

# -------------------------
# Create probabilities
# -------------------------

df["prob"] = df["trips"] / df["trips"].sum()

# -------------------------
# Group size model
# -------------------------

def sample_group_size(trips):
    if trips < 100:
        probs = [0.7, 0.2, 0.08, 0.02]
    elif trips < 500:
        probs = [0.5, 0.3, 0.15, 0.05]
    else:
        probs = [0.3, 0.4, 0.2, 0.1]

    return np.random.choice([1, 2, 3, 4], p=probs)

# -------------------------
# Settings
# -------------------------

DAYS = 5
GROUPS_PER_DAY = 60

START_TIME = 8 * 60
END_TIME   = 22 * 60

np.random.seed(42)

rows = []

# -------------------------
# Generate demand
# -------------------------

for day in range(1, DAYS + 1):

    # --- Demand-driven sampling ---
    sampled_indices = np.random.choice(
        df.index,
        size=GROUPS_PER_DAY,
        p=df["prob"],
        replace=True
    )

    for idx in sampled_indices:
        row = df.loc[idx]

        passengers = sample_group_size(row["trips"])
        time = np.random.randint(START_TIME, END_TIME)
        stopover = np.random.choice([0, 1])

        rows.append({
            "origin": kommune_to_vertiport[row["FraKommune"]],
            "destination": kommune_to_vertiport[row["TilKommune"]],
            "time": time,
            "number_of_passengers": passengers,
            "stopover_allowed": stopover,
            "day": day
        })

    # --- Exploration: ensure all vertiports are used ---
    used = set(
        [r["origin"] for r in rows if r["day"] == day] +
        [r["destination"] for r in rows if r["day"] == day]
    )

    missing = set(all_vertiports) - used

    for vp in missing:

        other = np.random.choice([v for v in all_vertiports if v != vp])

        if np.random.rand() < 0.5:
            origin, destination = vp, other
        else:
            origin, destination = other, vp

        time = np.random.randint(START_TIME, END_TIME)

        # small groups for random demand
        passengers = np.random.choice([1, 2], p=[0.8, 0.2])

        stopover = np.random.choice([0, 1])

        rows.append({
            "origin": origin,
            "destination": destination,
            "time": time,
            "number_of_passengers": passengers,
            "stopover_allowed": stopover,
            "day": day
        })

# -------------------------
# Create dataframe
# -------------------------

df_out = pd.DataFrame(rows)

# -------------------------
# Sort chronologically
# -------------------------

df_out = df_out.sort_values(["day", "time"]).reset_index(drop=True)

df_out["group"] = df_out.index + 1

df_out = df_out[
    [
        "group",
        "origin",
        "destination",
        "time",
        "number_of_passengers",
        "stopover_allowed",
        "day"
    ]
]

# -------------------------
# Save
# -------------------------

output_file = BASE_DIR / "synthetic_demand_groups.xlsx"
df_out.to_excel(output_file, index=False)

print(f"Saved: {output_file}")
print(df_out.head())