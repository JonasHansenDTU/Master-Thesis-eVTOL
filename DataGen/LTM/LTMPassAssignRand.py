import pandas as pd
import numpy as np

# -----------------------------
# SETTINGS
# -----------------------------
NUM_DAYS = 3
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

print(f"Saved to: {output_path}")
