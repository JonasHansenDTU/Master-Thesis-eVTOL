import pandas as pd
import matplotlib.pyplot as plt
import io

# 1. Load data
raw_data = """År    Måned    Middel minimum    Absolut Minimum    Døgn med frost
2026    Maj    7.1    -3.3    0.4
2026    April    3    -4.1    3.1
2026    Marts    1.9    -4.5    4.5
2026    Februar    -2.6    -16.8    21.2
2026    Januar    2.7    -19.4    24.5
2025    December    3.4    -7.5    6.3
2025    November    4.4    -9.1    8.3
2025    Oktober    7.2    -3.7    1
2025    September    10.8    -0.3    0
2025    August    12    4    0
2025    Juli    14.1    7.1    0
2025    Juni    10.6    1.2    0
2025    Maj    6.6    -1.9    0.2
2025    April    4.6    -4.6    2.6
2025    Marts    1.6    -8.7    8
2025    Februar    -0.7    -14    14.3
2025    Januar    0.4    -9.8    11.1
2024    December    3.3    -4.2    3.8
2024    November    3.3    -12    5.6
2024    Oktober    7.6    -1.7    0.1
2024    September    6.4    0.6    0
2024    August    12.7    4.6    0
2024    Juli    12.1    5.1    0
2024    Juni    10.4    1.6    0
2024    Maj    9.9    1.1    0
2024    April    4.2    -3.8    3.2
2024    Marts    2.8    -5.6    4.3
2024    Februar    2.4    -5.6    4.3
2024    Januar    -1.6    -16.4    16
2023    December    0.9    -16.5    10.3
2023    November    2    -15    9.4
2023    Oktober    7.1    -2.5    0.3
2023    September    12.6    1.8    0
2023    August    12.7    4.4    0
2023    Juli    11.9    5.5    0
2023    Juni    10.7    1.4    0
2023    Maj    6.4    -2.6    0.3
2023    April    2.9    -6.1    4.8
2023    Marts    0.2    -15.3    15.3
2023    Februar    1.1    -6.5    9.3
2023    Januar    1.7    -8.6    9
2022    December    -1.1    -18    16
2022    November    5.7    -6.2    3.6
2022    Oktober    8.9    -2.5    0.2
2022    September    9.7    -1.1    0
2022    August    13    4.3    0
2022    Juli    12    4    0
2022    Juni    10.8    2.5    0
2022    Maj    7.2    -2.3    0.3
2022    April    1.8    -8    6.6
2022    Marts    -0.7    -8.7    17.8
2022    Februar    1.4    -7.1    5.9
2022    Januar    1.9    -8.3    7.3
2021    December    -0.3    -17.5    15.7
2021    November    4.2    -5.2    3.6
2021    Oktober    7.5    -1.7    0.2
2021    September    10.9    2.7    0
2021    August    11.8    3.3    0
2021    Juli    14.1    6.6    0
2021    Juni    11.3    1.9    0
2021    Maj    6.5    -1.4    0.1
2021    April    1.1    -6.8    9.2
2021    Marts    0.2    -11.3    14.3
2021    Februar    -2.7    -20.7    16.6
2021    Januar    -1.2    -11.9    21.5
2020    December    2.6    -8.2    2.7
2020    November    5.3    -4.8    3.1
2020    Oktober    7.6    -2.6    0.3
2020    September    9.9    -1.3    0.1
2020    August    13.2    3.8    0
2020    Juli    19.8    1.8    0
2020    Juni    11.4    1.8    0
2020    Maj    5.6    -3.3    0.7
2020    April    2.8    -5.4    4.4
2020    Marts    0.9    -7.8    12.4
2020    Februar    2.3    -6.8    3.7
2020    Januar    3.3    -4.2    2.4"""

column_names = ['År', 'Måned', 'Middel minimum', 'Absolut Minimum', 'Døgn med frost']
df = pd.read_csv(io.StringIO(raw_data), sep='\s+', skiprows=1, names=column_names)

df = df.iloc[::-1].reset_index(drop=True)

months_en = {
    'Januar': 'Jan', 'Februar': 'Feb', 'Marts': 'Mar', 'April': 'Apr',
    'Maj': 'May', 'Juni': 'Jun', 'Juli': 'Jul', 'August': 'Aug',
    'September': 'Sep', 'Oktober': 'Oct', 'November': 'Nov', 'December': 'Dec'
}
df['Måned_en'] = df['Måned'].map(months_en)
df['Tid'] = df['Måned_en'] + " " + df['År'].astype(str)

# 2. Setup Plot
fig, ax1 = plt.subplots(figsize=(14, 6.5))
plt.rcParams['font.family'] = 'sans-serif'
ax1.grid(True, axis='y', linestyle='--', alpha=0.5, zorder=0)

color_bars = '#b0c4de'
color_line = '#8b0000'

bars = ax1.bar(df['Tid'], df['Døgn med frost'], color=color_bars, alpha=0.85,
               zorder=3, label='Frost Days')
ax1.set_ylabel('Number of Frost Days per Month', color='#4682b4',
               fontsize=12, fontweight='bold')
ax1.tick_params(axis='y', labelcolor='#4682b4', labelsize=11)
ax1.set_ylim(0, 31)

ax2 = ax1.twinx()
line, = ax2.plot(df['Tid'], df['Absolut Minimum'], color=color_line,
                 linewidth=2.5, marker='o', markersize=5, zorder=4,
                 label='Absolute Minimum Temp.')
ax2.set_ylabel('Temperature (°C)', color=color_line, fontsize=12, fontweight='bold')
ax2.tick_params(axis='y', labelcolor=color_line, labelsize=11)
ax2.set_ylim(-25, 10)

ax2.axhline(0, color='gray', linestyle='-', linewidth=1, alpha=0.5)
ax2.axhline(-8, color='darkorange', linestyle=':', linewidth=1.8, alpha=0.9,
            label='Very Cold Threshold (-8°C)')

# X-axis: show only January each year
ax1.set_xticks(range(len(df['Tid'])))
ax1.set_xticklabels(df['Tid'], rotation=45, ha='right', fontsize=10.5)

for index, label in enumerate(ax1.xaxis.get_ticklabels()):
    if df['Måned_en'].iloc[index] != 'Jan':
        label.set_visible(False)

lines, labels = ax1.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax1.legend(lines + lines2, labels + labels2, loc='upper left', frameon=True,
           facecolor='white', framealpha=0.9, fontsize=11)

plt.title('DMI Climate Data Analysis (2020-2026): Frost Days vs. Extreme Cold Events',
          fontsize=20, fontweight='bold', pad=15)
fig.tight_layout()

plt.savefig('dmi_temperature_analytics.png', dpi=300)
plt.show()