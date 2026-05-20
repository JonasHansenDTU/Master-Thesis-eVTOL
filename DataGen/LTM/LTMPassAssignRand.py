import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# SETTINGS
# -----------------------------
GROUPS_PER_DAY = 100

# Operating hours (minutes from midnight)
START_DAY = 0   # 07:00
END_DAY = 780    # 20:00

# -----------------------------
# LOAD DATA
# -----------------------------
BASE_PATH = "DataGen/LTM/"

# Load OD data
df = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"
)

# Load vertiport mapping
vertiports = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx",
    sheet_name="VertiportIDGiant"
)

# -----------------------------
# PREPROCESSING
# -----------------------------

# Remove self-loops
df = df[df["FraKommune"] != df["TilKommune"]]

# Use only air travel demand
df["demand"] = df["AntalErhvervsture_fly"]

# Remove OD pairs with zero demand
df = df[df["demand"] > 0]

# -----------------------------
# MAP MUNICIPALITIES TO VERTIPORTS
# -----------------------------

kommune_to_vertiport = dict(
    zip(vertiports["Kommuneid"], vertiports["id"])
)

df["origin_id"] = df["FraKommune"].map(kommune_to_vertiport)
df["dest_id"] = df["TilKommune"].map(kommune_to_vertiport)

# Remove rows without valid vertiports
df = df.dropna(subset=["origin_id", "dest_id"])

# Convert IDs to integers
df["origin_id"] = df["origin_id"].astype(int)
df["dest_id"] = df["dest_id"].astype(int)

# -----------------------------
# PROBABILITY DISTRIBUTION
# -----------------------------

# Higher observed demand -> higher sampling probability
probabilities = df["demand"] / df["demand"].sum()

# -----------------------------
# TIME SAMPLING
# -----------------------------

def sample_departure_time():
    """
    Creates realistic demand throughout the day:

    - Strong morning peak
    - Very low demand around lunch
    - Broader afternoon demand
    """

    r = np.random.rand()

    # -----------------------------
    # MORNING PEAK (45%)
    # -----------------------------
    if r < 0.45:

        # Around 09:00
        return int(np.random.normal(
            loc=120,
            scale=50
        ))

    # -----------------------------
    # AFTERNOON DEMAND (55%)
    # -----------------------------
    else:

        # Around 16:30
        return int(np.random.normal(
            loc=570,
            scale=110
        ))


def sample_departure_time_valid():
    """
    Keeps generating times until
    they fall within operating hours
    """

    while True:

        t = sample_departure_time()

        if START_DAY <= t <= END_DAY:
            return t

# -----------------------------
# PASSENGER GROUP SIZE
# -----------------------------

def sample_group_size():

    # Fully random group size
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

    # Extra safety
    if origin == destination:
        continue

    # Generate attributes
    departure_time = sample_departure_time_valid()

    passengers = sample_group_size()

    stopover = np.random.choice([0, 1])

    # Save generated trip
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

    # Random number of additional trips
    num_extra_trips = np.random.randint(1, 6)

    for _ in range(num_extra_trips):

        # Select random counterpart vertiport
        other = np.random.choice(
            list(all_vertiports - {vp})
        )

        # Randomize trip direction
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

# Sort chronologically
df_out = df_out.sort_values(
    by=["time"]
).reset_index(drop=True)

# Reassign IDs after sorting
df_out["group"] = range(1, len(df_out) + 1)

# -----------------------------
# SAVE OUTPUT
# -----------------------------

output_path = "inputData/LTM_demand.xlsx"

df_out.to_excel(output_path, index=False)

print(f"Saved to: {output_path}")

# -----------------------------
# VISUALIZE DEMAND
# -----------------------------

hours = df_out["time"] / 60

plt.figure(figsize=(10, 5))

plt.hist(hours, bins=20)

plt.xlabel("Time of Day")
plt.ylabel("Number of Passenger Groups")
plt.title("Synthetic eVTOL Demand Throughout the Day")

# Visual reference for lunch period
plt.axvline(12, linestyle="--")
plt.axvline(14, linestyle="--")

plt.show()