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
<<<<<<< HEAD
<<<<<<< HEAD
# Load OD data
=======
# Load data
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
# Load OD data
>>>>>>> d5dc66b (LTMDataGen)
# -------------------------

df = pd.read_excel(file_path, sheet_name="LTM")

<<<<<<< HEAD
<<<<<<< HEAD
MODE = "fly"
=======
MODE = "bil"
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
MODE = "fly"
>>>>>>> d5dc66b (LTMDataGen)
col = f"AntalErhvervsture_{MODE}"

df = df[["FraKommune", "TilKommune", col]].copy()
df.rename(columns={col: "trips"}, inplace=True)

<<<<<<< HEAD
<<<<<<< HEAD
# remove zero flows
df = df[df["trips"] > 0]

# remove self trips
df = df[df["FraKommune"] != df["TilKommune"]]

# -------------------------
# Load vertiports
# -------------------------

df_verti = pd.read_excel(file_path, sheet_name="VertiportID")

# mapping kommune → vertiport id
kommune_to_vertiport = dict(zip(df_verti["Kommuneid"], df_verti["id"]))

# keep only routes with vertiports
df = df[
    df["FraKommune"].isin(kommune_to_vertiport) &
    df["TilKommune"].isin(kommune_to_vertiport)
]

print(f"Ruter efter vertiport filter: {len(df)}")

# -------------------------
# Create probabilities
=======
# fjern null flows
=======
# remove zero flows
>>>>>>> d5dc66b (LTMDataGen)
df = df[df["trips"] > 0]

# remove self trips
df = df[df["FraKommune"] != df["TilKommune"]]

# -------------------------
<<<<<<< HEAD
# Create probabilities (VERY IMPORTANT)
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
# Load vertiports
# -------------------------

df_verti = pd.read_excel(file_path, sheet_name="VertiportID")

# mapping kommune → vertiport id
kommune_to_vertiport = dict(zip(df_verti["Kommuneid"], df_verti["id"]))

# keep only routes with vertiports
df = df[
    df["FraKommune"].isin(kommune_to_vertiport) &
    df["TilKommune"].isin(kommune_to_vertiport)
]

print(f"Ruter efter vertiport filter: {len(df)}")

# -------------------------
# Create probabilities
>>>>>>> d5dc66b (LTMDataGen)
# -------------------------

df["prob"] = df["trips"] / df["trips"].sum()

# -------------------------
<<<<<<< HEAD
<<<<<<< HEAD
# Group size model
=======
# Group size model (more group-heavy)
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
# Group size model
>>>>>>> d5dc66b (LTMDataGen)
# -------------------------

def sample_group_size(trips):
    if trips < 100:
        probs = [0.7, 0.2, 0.08, 0.02]
    elif trips < 500:
        probs = [0.5, 0.3, 0.15, 0.05]
    else:
        probs = [0.3, 0.4, 0.2, 0.1]

<<<<<<< HEAD
<<<<<<< HEAD
    return np.random.choice([1, 2, 3, 4], p=probs)
=======
    return np.random.choice([1,2,3,4], p=probs)
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
    return np.random.choice([1, 2, 3, 4], p=probs)
>>>>>>> d5dc66b (LTMDataGen)

# -------------------------
# Settings
# -------------------------

DAYS = 5
GROUPS_PER_DAY = 60

<<<<<<< HEAD
<<<<<<< HEAD
START_TIME = 8 * 60   # 08:00
END_TIME   = 22 * 60  # 22:00
=======
START_TIME = 8 * 60
END_TIME   = 22 * 60
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
START_TIME = 8 * 60   # 08:00
END_TIME   = 22 * 60  # 22:00
>>>>>>> d5dc66b (LTMDataGen)

np.random.seed(42)

rows = []

# -------------------------
<<<<<<< HEAD
<<<<<<< HEAD
# Generate synthetic demand
# -------------------------

for day in range(1, DAYS + 1):
=======
# Generate groups (KEY CHANGE)
# -------------------------

for day in range(1, DAYS+1):
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
# Generate synthetic demand
# -------------------------

for day in range(1, DAYS + 1):
>>>>>>> d5dc66b (LTMDataGen)

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
<<<<<<< HEAD
<<<<<<< HEAD
        stopover = np.random.choice([0, 1])

        rows.append({
            "origin": kommune_to_vertiport[row["FraKommune"]],
            "destination": kommune_to_vertiport[row["TilKommune"]],
=======
        stopover = np.random.choice([0,1])

        rows.append({
            "origin": row["FraKommune"],
            "destination": row["TilKommune"],
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
        stopover = np.random.choice([0, 1])

        rows.append({
            "origin": kommune_to_vertiport[row["FraKommune"]],
            "destination": kommune_to_vertiport[row["TilKommune"]],
>>>>>>> d5dc66b (LTMDataGen)
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
<<<<<<< HEAD
<<<<<<< HEAD
# Sort chronologically
=======
# Sort correctly
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
# Sort chronologically
>>>>>>> d5dc66b (LTMDataGen)
# -------------------------

df_out = df_out.sort_values(["day", "time"]).reset_index(drop=True)

<<<<<<< HEAD
<<<<<<< HEAD
# assign group ID AFTER sorting
df_out["group"] = df_out.index + 1

# reorder columns
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
=======
# assign group id AFTER sorting
df_out["group"] = df_out.index + 1

# reorder
df_out = df_out[[
    "group",
    "origin",
    "destination",
    "time",
    "number_of_passengers",
    "stopover_allowed",
    "day"
]]
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
# assign group ID AFTER sorting
df_out["group"] = df_out.index + 1

# reorder columns
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
>>>>>>> d5dc66b (LTMDataGen)

# -------------------------
# Save
# -------------------------

output_file = BASE_DIR / "synthetic_demand_groups.xlsx"
df_out.to_excel(output_file, index=False)

<<<<<<< HEAD
<<<<<<< HEAD
print(f"\n✅ Saved: {output_file}")
=======
print(f"✅ Saved: {output_file}")
>>>>>>> f35fef3 (LTM Data to PassengerGroups)
=======
print(f"\n✅ Saved: {output_file}")
>>>>>>> d5dc66b (LTMDataGen)
print(df_out.head())