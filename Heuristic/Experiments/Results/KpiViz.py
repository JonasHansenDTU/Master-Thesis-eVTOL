import pandas as pd
import matplotlib.pyplot as plt
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))     
SAVE_PATH  = os.path.join(SCRIPT_DIR, 'KPIsViz.png')
SUMMARY_CSV = os.path.join(SCRIPT_DIR, 'Scenarios', 'scenario_summary.csv')

# 1. Load computed averages and stds from ScenarioSummary.py output
raw = pd.read_csv(SUMMARY_CSV)

SCENARIO_ORDER = ['Base', 's1', 's2', 's3', 's4', 's5', 's6', 's7', 's8', 's9', 's10']
raw['_order'] = raw['scenario'].map({s: i for i, s in enumerate(SCENARIO_ORDER)})
raw = raw.sort_values('_order').drop(columns='_order').reset_index(drop=True)

# Pretty x-axis labels
raw['label'] = raw['scenario'].replace({'Base': 'Base'}).apply(
    lambda s: s if s == 'Base' else 'S' + s[1:]
)

# 2. Metric definitions: (mean_col, std_col, y-axis title)
metrics = [
    ('best_obj_mean',                'best_obj_std',                'Objective Value'),
    ('best_profit_mean',             'best_profit_std',             'Profit (DKK)'),
    ('aircraft_utilization_mean',    'aircraft_utilization_std',    'Aircraft Utilization'),
    ('deadhead_time_mean',           'deadhead_time_std',           'Deadhead Time (min)'),
    ('passenger_demand_served_mean', 'passenger_demand_served_std', '# Passengers Served'),
]

# 3. Plot — 2×3 grid, 5 panels used
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle(
    f'Sensitivity Analysis: KPI Comparisons Across {len(raw)} Scenarios (Mean ± SD)',
    fontsize=16, fontweight='bold', y=0.96
)
axes = axes.flatten()

for ax, (mean_col, sd_col, title) in zip(axes, metrics):
    ax.errorbar(
        raw['label'], raw[mean_col], yerr=raw[sd_col],
        fmt='o', color='#C95F5F', ecolor='#AAAAAA',
        elinewidth=2, capsize=4, markersize=7
    )

    if raw[mean_col].min() - raw[sd_col].max() < 0:
        ax.axhline(0, color='black', linestyle='--', alpha=0.3, lw=1.2)

    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_axisbelow(True)
    ax.tick_params(axis='x', rotation=30)

# Remove the empty 6th panel
fig.delaxes(axes[5])

# 4. Export
plt.tight_layout(rect=[0, 0.03, 1, 0.92])
plt.savefig(SAVE_PATH, dpi=300, bbox_inches='tight')
plt.show()