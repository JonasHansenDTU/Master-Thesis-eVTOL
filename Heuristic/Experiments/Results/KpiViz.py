import pandas as pd
import matplotlib.pyplot as plt
import os


SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(SCRIPT_DIR, 'KPIsViz.png')

# 1. Complete dataset for all 10 scenarios (decimal commas converted to periods)
data = {
    'Scenario': ['Sc. 1', 'Sc. 2', 'Sc. 3', 'Sc. 4', 'Sc. 5', 'Sc. 6', 'Sc. 7', 'Sc. 8', 'Sc. 9', 'Sc. 10'],
    'Obj. Value': [-3695, 6698, -6334, -6769, 8648, 651, -9614, 10261, 1729, 9306],
    'Obj. Value SD': [7243, 8597, 5836, 5156, 9553, 8676, 3196, 9149, 3718, 9339],
    'Profit': [12876, 23352, 10710, 9882, 25205, 17552, 7354, 26402, 18245, 25297],
    'Profit SD': [6987, 8723, 5554, 4757, 8975, 8470, 3149, 9005, 3744, 9349],
    'Aircraft Uti.': [0.11, 0.12, 0.20, 0.14, 0.09, 0.17, 0.17, 0.13, 0.18, 0.14],
    'Aircraft Uti. SD': [0.03, 0.01, 0.07, 0.02, 0.01, 0.05, 0.06, 0.05, 0.05, 0.05],
    'Deadhead': [41.6, 115.0, 17.5, 44.3, 132.1, 66.2, 11.7, 93.1, 126.4, 72.4],
    'Deadhead SD': [36.7, 70.7, 20.7, 43.4, 86.5, 48.1, 30.7, 47.7, 93.5, 51.0],
    'Passengers': [5.1, 7.4, 4.0, 5.0, 7.3, 5.2, 3.1, 7.4, 6.6, 7.6],
    'Passengers SD': [2.3, 1.8, 1.8, 2.7, 2.1, 2.6, 1.8, 2.1, 2.4, 2.2]
}

df = pd.DataFrame(data)

# 2. Set up a clean grid layout for the 5 metrics (2 rows, 3 columns)
fig, axes = plt.subplots(2, 3, figsize=(18, 11))
fig.suptitle('Sensitivity Analysis: Metric Comparisons Across 10 Scenarios (Mean ± SD)', fontsize=16, fontweight='bold', y=0.96)

# Define columns, error bounds, and titles
metrics = [
    ('Obj. Value', 'Obj. Value SD', 'Objective Value'),
    ('Profit', 'Profit SD', 'Profit ($)'),
    ('Aircraft Uti.', 'Aircraft Uti. SD', 'Aircraft Utilization'),
    ('Deadhead', 'Deadhead SD', 'Deadhead Miles/Time'),
    ('Passengers', 'Passengers SD', '# of Passengers Served')
]

axes = axes.flatten()

# 3. Generate plots dynamically loop through metrics
for i, (mean_col, sd_col, title) in enumerate(metrics):
    ax = axes[i]
    
    # Blue dot markers with muted grey error bars for optimal legibility
    ax.errorbar(df['Scenario'], df[mean_col], yerr=df[sd_col], 
                fmt='o', color='#1f77b4', ecolor='#95a5a6', 
                elinewidth=2, capsize=4, markersize=7)
    
    # Add a horizontal baseline indicator at 0 if the data has negative entries
    if df[mean_col].min() - df[sd_col].max() < 0:
        ax.axhline(0, color='black', linestyle='--', alpha=0.3, lw=1.2)
        
    ax.set_title(title, fontsize=12, fontweight='bold', pad=10)
    ax.grid(True, linestyle=':', alpha=0.6)
    ax.set_axisbelow(True)
    
    # Rotate the scenario text markers slightly to prevent 10 columns from clustering
    ax.tick_params(axis='x', rotation=30)

# 4. Clean up the grid by dropping the empty 6th panel
fig.delaxes(axes[5])

# 5. Handle formatting boundary spaces and export to file
plt.tight_layout(rect=[0, 0.03, 1, 0.92])
plt.savefig(SAVE_PATH, dpi=300, bbox_inches='tight')
plt.show()