import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# SETTINGS
# -----------------------------
GROUPS_PER_DAY = 100

# Operating windows (minutes from midnight)
START_MORNING = 420    # 07:00
END_MORNING = 720      # 12:00

START_AFTERNOON = 840  # 14:00
END_DAY = 1200         # 20:00

# -----------------------------
# LOAD DATA
# -----------------------------
BASE_PATH = "DataGen/LTM/"

# Load OD data
df = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"
)

# Load vertiport sheet
vertiports = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx",
    sheet_name="VertiportID"
)

# -----------------------------
# PREPROCESS
# -----------------------------

# Remove trips from municipality to itself
df = df[df["FraKommune"] != df["TilKommune"]]

# Only use air travel demand
df["demand"] = df["AntalErhvervsture_fly"]

# Remove OD pairs with zero demand
df = df[df["demand"] > 0]

# -----------------------------
# MAP MUNICIPALITY -> VERTIPORT
# -----------------------------

# Create dictionary for mapping
kommune_to_vertiport = dict(
    zip(vertiports["Kommuneid"], vertiports["id"])
)

# Map municipality IDs to vertiport IDs
df["origin_id"] = df["FraKommune"].map(kommune_to_vertiport)
df["dest_id"] = df["TilKommune"].map(kommune_to_vertiport)

# Remove rows without valid vertiports
df = df.dropna(subset=["origin_id", "dest_id"])

# Convert IDs to integers
df["origin_id"] = df["origin_id"].astype(int)
df["dest_id"] = df["dest_id"].astype(int)

# -----------------------------
# DEMAND PROBABILITIES
# -----------------------------

# Higher demand = higher probability
probabilities = df["demand"] / df["demand"].sum()

# -----------------------------
# TIME SAMPLING
# -----------------------------

def sample_departure_time():
    """
    Generates realistic departure times.

    Morning:
    - Shorter operating window
    - More concentrated demand

    Afternoon:
    - Longer operating window
    - More spread-out demand
    """

    r = np.random.rand()

    # 45% morning departures
    if r < 0.45:

        # Peak around 09:00
        return int(np.random.normal(loc=540, scale=45))

    # 55% afternoon departures
    else:

        # Peak around 16:30
        return int(np.random.normal(loc=990, scale=140))


def sample_departure_time_valid():
    """
    Continues generating departure times
    until they fall within valid operating hours
    """

    while True:

        t = sample_departure_time()

        if (
            START_MORNING <= t <= END_MORNING
        ) or (
            START_AFTERNOON <= t <= END_DAY
        ):
            return t

# -----------------------------
# GROUP SIZE
# -----------------------------

def sample_group_size():
    """
    Random passenger group size
    between 1 and 4 passengers
    """

    return np.random.choice([1, 2, 3, 4])

# -----------------------------
# GENERATE DEMAND
# -----------------------------

rows = []
group_id = 1

# Sample OD pairs based on probabilities
sampled_indices = np.random.choice(
    df.index,
    size=GROUPS_PER_DAY,
    p=probabilities
)

used_vertiports = set()

for idx in sampled_indices:

    row = df.loc[idx]

    origin = row["origin_id"]
    destination = row["dest_id"]

    # Extra safety check
    if origin == destination:
        continue

    # Generate trip attributes
    departure_time = sample_departure_time_valid()

    passengers = sample_group_size()

    stopover = np.random.choice([0, 1])

    # Save trip
    rows.append({
        "group": group_id,
        "origin": origin,
        "destination": destination,
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

    # Random number of extra trips
    num_extra_trips = np.random.randint(1, 6)

    for _ in range(num_extra_trips):

        # Random counterpart vertiport
        other = np.random.choice(list(all_vertiports - {vp}))

        # Random trip direction
        if np.random.rand() < 0.5:
            origin = vp
            destination = other
        else:
            origin = other
            destination = vp

        rows.append({
            "group": group_id,
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

# Sort chronologically by departure time
df_out = df_out.sort_values(
    by=["time"]
).reset_index(drop=True)

# Reassign group IDs after sorting
df_out["group"] = range(1, len(df_out) + 1)

# -----------------------------
# SAVE EXCEL FILE
# -----------------------------

output_path = "inputData/LTM_demand.xlsx"

df_out.to_excel(output_path, index=False)

print(f"Saved to: {output_path}")

# -----------------------------
# DEMAND VISUALIZATION
# -----------------------------

# Convert minutes to hours
hours = df_out["time"] / 60

plt.figure(figsize=(10, 5))

plt.hist(hours, bins=20)

plt.xlabel("Time of Day")
plt.ylabel("Number of Passenger Groups")
plt.title("Synthetic eVTOL Demand Throughout the Day")

# Visualize break period
plt.axvline(12, linestyle="--")
plt.axvline(14, linestyle="--")

plt.show()