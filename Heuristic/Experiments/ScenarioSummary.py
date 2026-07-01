import pandas as pd
from pathlib import Path

KPI_COLS = [
    "best_obj",
    "best_profit",
    "aircraft_utilization",
    "deadhead_time",
    "passenger_demand_served",
    "average_passenger_waiting_time",
]

SCENARIO_ORDER = ["Base", "s1", "s2", "s3", "s4", "s5", "s6", "s7", "s8", "s9", "s10"]

results_dir = Path(__file__).parent / "Results" / "Scenarios"
df = pd.read_csv(results_dir / "all_scenarios_kpis.csv")

agg = (
    df.groupby("scenario")[KPI_COLS]
    .agg(["mean", "std"])
)

# Flatten multi-level columns to e.g. "best_obj_mean", "best_obj_std"
agg.columns = [f"{col}_{stat}" for col, stat in agg.columns]
agg = agg.reset_index()

# Reorder rows to match scenario order
agg["_order"] = agg["scenario"].map(
    {s: i for i, s in enumerate(SCENARIO_ORDER)}
)
agg = agg.sort_values("_order").drop(columns="_order").reset_index(drop=True)

out_path = results_dir / "scenario_summary.csv"
agg.to_csv(out_path, index=False)
print(f"Saved summary to: {out_path}")
print()
print(agg.to_string(index=False))
