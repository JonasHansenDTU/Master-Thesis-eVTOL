<<<<<<< HEAD
<<<<<<< HEAD
<<<<<<< HEAD
# %%
=======
<<<<<<< HEAD
>>>>>>> 20220b3 (LTMDataGen)
=======
<<<<<<< HEAD
>>>>>>> 8c31caa (Updated files)
=======
<<<<<<< HEAD
>>>>>>> 2b422b6 (Fix)
import pandas as pd
import numpy as np
from pathlib import Path

# -------------------------
# Paths
# -------------------------

<<<<<<< HEAD
BASE_DIR = Path("/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM")
file_path = BASE_DIR / "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"

# -------------------------
=======
# Operating windows (minutes from midnight)
START_MORNING = 420    # 07:00
END_MORNING = 720      # 12:00

START_AFTERNOON = 840  # 14:00
END_DAY = 1200         # 20:00
=======
LTMPassAssignRand.py
# %%
=======
>>>>>>> 04c85b4 (Fix)
import pandas as pd
import numpy as np

<<<<<<< HEAD
# -------------------------
# Paths
# -------------------------

BASE_DIR = Path("/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM")
file_path = BASE_DIR / "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"
>>>>>>> ea106b8 (Updated files)

# -------------------------
# Load OD data
# -------------------------

<<<<<<< HEAD
>>>>>>> 8c31caa (Updated files)
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
<<<<<<< HEAD
        probs = [0.3, 0.4, 0.2, 0.1]

    return np.random.choice([1, 2, 3, 4], p=probs)

# -------------------------
# Settings
# -------------------------
=======
        # Peak around 16:30
        return int(np.random.normal(loc=990, scale=140))
=======
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
>>>>>>> ea106b8 (Updated files)

    return np.random.choice([1, 2, 3, 4], p=probs)

<<<<<<< HEAD
def sample_departure_time_valid():
    """
    Keeps generating times until
    they fall within operating windows
    """
>>>>>>> 8c31caa (Updated files)

DAYS = 5
GROUPS_PER_DAY = 60

START_TIME = 8 * 60
END_TIME   = 22 * 60

np.random.seed(42)

=======
# -------------------------
# Settings
# -------------------------

DAYS = 5
GROUPS_PER_DAY = 60

START_TIME = 8 * 60
END_TIME   = 22 * 60

np.random.seed(42)

>>>>>>> ea106b8 (Updated files)
rows = []

# -------------------------
# Generate demand
# -------------------------

<<<<<<< HEAD
for day in range(1, DAYS + 1):

    # --- Demand-driven sampling ---
=======
<<<<<<< HEAD
    # Sample OD pairs based on probabilities
>>>>>>> 8c31caa (Updated files)
    sampled_indices = np.random.choice(
        df.index,
        size=GROUPS_PER_DAY,
        p=df["prob"],
        replace=True
    )
=======
for day in range(1, DAYS + 1):
>>>>>>> ea106b8 (Updated files)

<<<<<<< HEAD
=======
    # --- Demand-driven sampling ---
    sampled_indices = np.random.choice(
        df.index,
        size=GROUPS_PER_DAY,
        p=df["prob"],
        replace=True
    )

>>>>>>> 8c31caa (Updated files)
    for idx in sampled_indices:

        row = df.loc[idx]

<<<<<<< HEAD
        origin = row["origin_id"]
        destination = row["dest_id"]

        # Extra safety check
        if origin == destination:
            continue

        # Generate attributes
        departure_time = sample_departure_time_valid()
        passengers = sample_group_size()
        stopover = np.random.choice([0, 1])

        # Save trip
=======
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

>>>>>>> ea106b8 (Updated files)
        rows.append({
            "origin": origin,
            "destination": destination,
<<<<<<< HEAD
            "time": departure_time,
            "number_of_passengers": passengers,
            "stopover_allowed": stopover
        })

        used_vertiports.add(origin)
        used_vertiports.add(destination)

        group_id += 1

    # -----------------------------
    # NETWORK COVERAGE ADJUSTMENT
    # -----------------------------

    all_vertiports = set(vertiports["id"])

    missing = all_vertiports - used_vertiports

    for vp in missing:

        # Random number of additional trips
        num_extra_trips = np.random.randint(1, 6)

        for _ in range(num_extra_trips):

            # Random other vertiport
            other = np.random.choice(list(all_vertiports - {vp}))

            # Random direction
            if np.random.rand() < 0.5:
                origin = vp
                destination = other
            else:
                origin = other
                destination = vp

            rows.append({
                "group": group_id,
                "day": day,
                "origin": origin,
                "destination": destination,
                "time": sample_departure_time_valid(),
                "number_of_passengers": sample_group_size(),
                "stopover_allowed": np.random.choice([0, 1])
            })

            group_id += 1

# -----------------------------
# FINAL DATAFRAME
# -----------------------------

df_out = pd.DataFrame(rows)

# Sort chronologically
df_out = df_out.sort_values(
    by=["day", "time"]
).reset_index(drop=True)

# Reassign IDs after sorting
df_out["group"] = range(1, len(df_out) + 1)

# -----------------------------
# SAVE EXCEL FILE
# -----------------------------

output_path = BASE_PATH + "synthetic_demand.xlsx"

df_out.to_excel(output_path, index=False)

print(f"Saved to: {output_path}")
<<<<<<< HEAD
=======
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
>>>>>>> d5dc66b (LTMDataGen)
=======

# -----------------------------
# DEMAND VISUALIZATION
# -----------------------------

# Convert minutes -> hours
hours = df_out["time"] / 60

plt.figure(figsize=(10,5))

plt.hist(hours, bins=20)

plt.xlabel("Time of Day")
plt.ylabel("Number of Passenger Groups")
plt.title("Synthetic eVTOL Demand Throughout the Day")

# Add vertical lines for operating windows
plt.axvline(12, linestyle="--")
plt.axvline(14, linestyle="--")

plt.show()
<<<<<<< HEAD
>>>>>>> 0cadae1 (More Liquid demand)
=======
=======
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
=======
# -----------------------------
# SETTINGS
# -----------------------------
NUM_DAYS = 3
>>>>>>> 04c85b4 (Fix)
GROUPS_PER_DAY = 60

START_MORNING = 420   # 07:00
END_MORNING   = 720   # 12:00
START_AFTERNOON = 840 # 14:00
END_DAY       = 1200  # 20:00

# -----------------------------
# LOAD DATA
# -----------------------------
BASE_PATH = "/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/DataGen/LTM/"

df = pd.read_excel(BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx")

vertiports = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx",
    sheet_name="VertiportID"
)

# -----------------------------
# PREPROCESS
# -----------------------------
df = df[df["FraKommune"] != df["TilKommune"]]

df["demand"] = df["AntalErhvervsture_fly"]
df = df[df["demand"] > 0]

# -----------------------------
# MAP kommune → vertiport
# -----------------------------
kommune_to_vertiport = dict(zip(vertiports["Kommuneid"], vertiports["id"]))

df["origin_id"] = df["FraKommune"].map(kommune_to_vertiport)
df["dest_id"]   = df["TilKommune"].map(kommune_to_vertiport)

df = df.dropna(subset=["origin_id", "dest_id"])

df["origin_id"] = df["origin_id"].astype(int)
df["dest_id"]   = df["dest_id"].astype(int)

# -----------------------------
# PROBABILITIES
# -----------------------------
probabilities = df["demand"] / df["demand"].sum()

# -----------------------------
# TIME SAMPLING
# -----------------------------
def sample_departure_time():
    r = np.random.rand()

    if r < 0.6:
        return int(np.random.normal(loc=600, scale=90))   # morning peak
    else:
        return int(np.random.normal(loc=1020, scale=120)) # afternoon peak


def sample_departure_time_valid():
    while True:
        t = sample_departure_time()
        if (START_MORNING <= t <= END_MORNING) or (START_AFTERNOON <= t <= END_DAY):
            return t

# -----------------------------
# GROUP SIZE (FULLY RANDOM)
# -----------------------------
def sample_group_size():
    return np.random.choice([1, 2, 3, 4])  # equal probability

# -----------------------------
# GENERATE DEMAND
# -----------------------------
rows = []
group_id = 1

for day in range(1, NUM_DAYS + 1):

    sampled_indices = np.random.choice(df.index, size=GROUPS_PER_DAY, p=probabilities)

    used_vertiports = set()

    for idx in sampled_indices:
        row = df.loc[idx]

        origin = row["origin_id"]
        dest   = row["dest_id"]

        if origin == dest:
            continue

        time = sample_departure_time_valid()
        passengers = sample_group_size()
        stopover = np.random.choice([0, 1])

        rows.append({
            "group": group_id,
            "day": day,
            "origin": origin,
            "destination": dest,
            "time": int(time),
            "number_of_passengers": int(passengers),
            "stopover_allowed": int(stopover)
        })

        used_vertiports.add(origin)
        used_vertiports.add(dest)

        group_id += 1

    # -----------------------------
    # ENSURE ALL VERTIPORTS USED
    # -----------------------------
    all_vertiports = set(vertiports["id"])
    missing = all_vertiports - used_vertiports

for vp in missing:
    # Random number of trips (1 to 5)
    num_extra_trips = np.random.randint(1, 6)

    for _ in range(num_extra_trips):

        other = np.random.choice(list(all_vertiports - {vp}))

        # Randomly decide direction (from or to vp)
        if np.random.rand() < 0.5:
            origin = vp
            destination = other
        else:
            origin = other
            destination = vp

        rows.append({
            "group": group_id,
            "day": day,
            "origin": origin,
            "destination": destination,
            "time": sample_departure_time_valid(),
            "number_of_passengers": np.random.choice([1, 2, 3, 4]),
            "stopover_allowed": np.random.choice([0, 1])
        })

        group_id += 1

# -----------------------------
# FINAL DATAFRAME
# -----------------------------
df_out = pd.DataFrame(rows)

df_out = df_out.sort_values(by=["day", "time"]).reset_index(drop=True)
df_out["group"] = range(1, len(df_out) + 1)

# -----------------------------
# SAVE
# -----------------------------
output_path = BASE_PATH + "synthetic_demand.xlsx"
df_out.to_excel(output_path, index=False)

<<<<<<< HEAD
=======
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

>>>>>>> ea106b8 (Updated files)
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
<<<<<<< HEAD
print(df_out.head())
>>>>>>> 81938cb (LTMDataGen)
<<<<<<< HEAD
>>>>>>> 20220b3 (LTMDataGen)
=======
=======
print(df_out.head())
>>>>>>> ea106b8 (Updated files)
<<<<<<< HEAD
>>>>>>> 8c31caa (Updated files)
=======
=======
print(f"Saved to: {output_path}")
>>>>>>> 04c85b4 (Fix)
>>>>>>> 2b422b6 (Fix)
