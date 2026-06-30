import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# SETTINGS
# =====================================================

GROUPS_PER_DAY = 50

# -----------------------------------------------------
# SIMULATION TIME
#
# IMPORTANT:
# 0 minutes = 07:00
# 780 minutes = 20:00
#
# This makes the simulation start at t = 0
# instead of t = 420.
# =====================================================

START_DAY = 0
END_DAY = 780

# =====================================================
# LOAD DATA
# =====================================================

BASE_PATH = "DataGen/LTM/"

# Load OD demand data
df = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx"
)

# Load vertiport mapping sheet
vertiports = pd.read_excel(
    BASE_PATH + "AntalErhvervstureMellemKommuner_GMM_Basis2025.xlsx",
    sheet_name="VertiportID"
)

# =====================================================
# PREPROCESSING
# =====================================================

# Remove trips where origin = destination
df = df[df["FraKommune"] != df["TilKommune"]]

# Use ONLY air travel demand
df["demand"] = df["AntalErhvervsture_fly"]

# Remove OD pairs with zero demand
df = df[df["demand"] > 0]

# =====================================================
# MAP MUNICIPALITIES -> VERTIPORTS
# =====================================================

# Create lookup dictionary
kommune_to_vertiport = dict(
    zip(vertiports["Kommuneid"], vertiports["id"])
)

# Map municipality codes to vertiport IDs
df["origin_id"] = df["FraKommune"].map(
    kommune_to_vertiport
)

df["dest_id"] = df["TilKommune"].map(
    kommune_to_vertiport
)

# Remove rows without valid vertiports
df = df.dropna(subset=["origin_id", "dest_id"])

# Convert IDs to integers
df["origin_id"] = df["origin_id"].astype(int)
df["dest_id"] = df["dest_id"].astype(int)

# =====================================================
# PROBABILITY DISTRIBUTION
# =====================================================

# Higher demand -> higher probability
probabilities = (
    df["demand"] / df["demand"].sum()
)

# =====================================================
# TIME SAMPLING
# =====================================================

def sample_departure_time():
    """
    Generates realistic departure times.

    Simulation time:
    0 = 07:00

    Structure:
    - Strong morning peak
    - Very low lunch demand
    - Broader afternoon demand
    """

    r = np.random.rand()

    # -------------------------------------------------
    # MORNING PEAK (45%)
    # -------------------------------------------------
    if r < 0.45:

        # Around 09:00
        # 09:00 -> 120 minutes after 07:00
        return int(
            np.random.normal(
                loc=120,
                scale=50
            )
        )

    # -------------------------------------------------
    # AFTERNOON DEMAND (55%)
    # -------------------------------------------------
    else:

        # Around 16:30
        # 16:30 -> 570 minutes after 07:00
        return int(
            np.random.normal(
                loc=570,
                scale=110
            )
        )


def sample_departure_time_valid():
    """
    Keeps generating times until
    they fall within operating hours
    """

    while True:

        t = sample_departure_time()

        if START_DAY <= t <= END_DAY:
            return t

# =====================================================
# PASSENGER GROUP SIZE
# =====================================================

def sample_group_size():

    # Fully random group size
    return np.random.choice([1, 2, 3, 4], p=[0.1, 0.2, 0.3, 0.4])

# =====================================================
# GENERATE DEMAND
# =====================================================

rows = []
group_id = 1

# Sample OD pairs
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

# =====================================================
# NETWORK COVERAGE ADJUSTMENT
# =====================================================

all_vertiports = set(
    vertiports["id"]
)

missing = all_vertiports - used_vertiports

for vp in missing:

    # Generate random number of extra trips
    num_extra_trips = np.random.randint(1, 6)

    for _ in range(num_extra_trips):

        # Random counterpart vertiport
        other = np.random.choice(
            list(all_vertiports - {vp})
        )

        # Randomize direction
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

# =====================================================
# FINAL DATAFRAME
# =====================================================

df_out = pd.DataFrame(rows)

# Sort chronologically
df_out = df_out.sort_values(
    by=["time"]
).reset_index(drop=True)

# Reassign group IDs after sorting
df_out["group"] = range(
    1,
    len(df_out) + 1
)

# =====================================================
# SAVE OUTPUT
# =====================================================

output_path = "inputData/LTM_demand.xlsx"
output_path = "inputData/LTM_demandBFinal.xlsx"

df_out.to_excel(
    output_path,
    index=False
)

print(f"Saved to: {output_path}")

# =====================================================
# DEMAND VISUALIZATION
# =====================================================

# Convert simulation time back to clock time
# because 0 = 07:00
# hours = (df_out["time"] + 420) / 60

# plt.figure(figsize=(10, 5))

# plt.hist(hours, bins=20)

# plt.xlabel("Time of Day")
# plt.ylabel("Number of Passenger Groups")

# plt.title(
#     "Synthetic eVTOL Demand Throughout the Day"
# )

# # Show proper clock hours
# plt.xticks(
#     [7,8,9,10,11,12,13,14,15,16,17,18,19,20]
# )

# # Visual reference for lunch hours
# plt.axvline(12, linestyle="--")
# plt.axvline(14, linestyle="--")

# plt.show()