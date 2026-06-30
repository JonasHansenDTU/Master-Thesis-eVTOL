import os
import re
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(SCRIPT_DIR, 'evtol_schedule_english.png')

# ==========================================
# THEME (matches the other thesis figures)
# ==========================================
LIGHT_RED = '#F9EAEA'
MID_RED   = '#C95F5F'
DARK_RED  = '#641818'
GRID_GRAY = '#AAAAAA'

plt.rcParams.update({"font.family": "serif"})

# One red per plane, all within the theme's red family but varied in
# lightness and hue. The order alternates light and dark shades so that
# two planes in adjacent rows never share a similar colour. Cycles if a
# panel has more planes than colours.
PLANE_COLORS = [
    '#F0B8B0',  # light salmon
    '#863232',  # dark brick
    '#E59A9A',  # dusty rose
    '#9A3B3B',  # deep brick
    '#D9806E',  # coral red
    '#B5484A',  # brick red
    '#CE6A6A',  # rose red
    '#A85C5C',  # muted mauve-red
    '#C0744F',  # terracotta
    '#C95F5F',  # theme mid red
]

# ==========================================
# MASTER DATA CONFIGURATION
# ==========================================
all_solutions = [
    {
        "title": "Base case (Max Profit: 25,432.53)",
        "flights": {
            'Plane 3':  [(153, 180, 'Group 21 (q=4)\n2->3'), (611, 638, 'Group 52 (q=4)\n3->2')],
            'Plane 5':  [(106, 119, 'Group 7 (q=1)\n3->5'),  (149, 162, 'Group 18 (q=3)\n5->3')],
            'Plane 10': [(216, 230, 'Group 27 (q=3)\n5->6'), (260, 274, 'Empty\n6->5')],
            'Plane 20': [(92, 117, 'Group 2 (q=4)\n10->4'),  (722, 752, 'Group 60 (q=4)\n4->1')]
        }
    },
    {
        "title": "Scenario 1 (Max Profit: 31,820.06)",
        "flights": {
            'Plane 1':  [(105, 115, 'Group 7 (q=3)\n1->10'), (622, 641, 'Empty\n10->3'), (656, 674, 'Empty\n3->1')],
            'Plane 5':  [(395, 411, 'Group 26 (q=3)\n3->4'), (477, 507, 'Group 29 (q=4)\n4->1'), (522, 540, 'Empty\n1->3')],
            'Plane 19': [(15, 25, 'Group 46 (q=3)\n10->3'), (79, 109, 'Group 5 (q=3)\n1->4'), (494, 524, 'Group 30 (q=4)\n4->1'), (560, 590, 'Group 39 (q=4)\n1->4'), (689, 719, 'Group 54 (q=4)\n4->1')]
        }
    },
    {
        "title": "Scenario 2 (Max Profit: 43,897.16)",
        "flights": {
            'Plane 1':  [(586, 616, 'Group 42 (q=3)\n1->4'), (646, 676, 'Group 49 (q=4)\n4->1')],
            'Plane 6':  [(0, 13, 'Empty\n3->5'), (607, 620, 'Group 43 (q=2)\n5->3')],
            'Plane 8':  [(145, 175, 'Group 20 (q=3)\n4->1'), (741, 771, 'Group 58 (q=4)\n1->4')],
            'Plane 10': [(609, 622, 'Empty\n5->3'), (652, 665, 'Empty\n3->5')],
            'Plane 15': [(42, 77, 'Group 2 (q=4)\n8->10'), (470, 505, 'Group 33 (q=2)\n10->8')],
            'Plane 18': [(122, 150, 'Group 17 (q=4)\n9->3'), (600, 628, 'Empty\n3->9')]
        }
    },
    {
        "title": "Scenario 3 (Max Profit: 21,064.12)",
        "flights": {
            'Plane 1':  [(30, 59, 'Group 3 (q=4)\n1->8'), (201, 236, 'Group 22 (q=4)\n8->10')],
            'Plane 3':  [(30, 46, 'Empty\n2->5')],
            'Plane 5':  [(513, 532, 'Group 35 (q=4)\n3->10')],
            'Plane 9':  [(520, 534, 'Group 36 (q=4)\n5->6')],
            'Plane 11': [(705, 731, 'Group 56 (q=4)\n6->3'), (761, 774, 'Empty\n3->5')],
            'Plane 20': [(687, 697, 'Group 55 (q=4)\n10->1')]
        }
    },
    {
        "title": "Scenario 4 (Max Profit: 18,041.75)",
        "flights": {
            'Plane 7':  [(42, 72, 'Group 3 (q=4)\n4->1'), (102, 132, 'Group 14 (q=4)\n1->4')],
            'Plane 8':  [(127, 157, 'Group 21 (q=3)\n4->1'), (358, 388, 'Group 30 (q=2)\n1->4')]
        }
    },
    {
        "title": "Scenario 5 (Max Profit: 38,737.11)",
        "flights": {
            'Plane 1':  [(752, 762, 'Group 62 (q=3)\n1->10')],
            'Plane 7':  [(71, 96, 'Group 4 (q=4)\n4->10'), (174, 193, 'Group 25 (q=3)\n10->3'), (208, 224, 'Empty\n3->4')],
            'Plane 9':  [(112, 125, 'Group 11 (q=4)\n5->3'), (628, 641, 'Group 47 (q=4)\n3->5')],
            'Plane 10': [(655, 668, 'Group 53 (q=4)\n5->3'), (758, 771, 'Empty\n3->5')],
            'Plane 12': [(638, 652, 'Group 49 (q=3)\n6->5'), (667, 681, 'Empty\n5->6')],
            'Plane 19': [(15, 25, 'Empty\n10->1'), (657, 687, 'Group 54 (q=4)\n1->4'), (720, 750, 'Group 58 (q=4)\n4->1')]
        }
    },
    {
        "title": "Scenario 6 (Max Profit: 46,153.42)",
        "flights": {
            'Plane 1':  [(62, 72, 'Group 5 (q=4)\n1->10'), (239, 249, 'Empty\n10->1'), (352, 362, 'Empty\n1->10')],
            'Plane 2':  [(181, 211, 'Group 32 (q=3)\n1->5'), (607, 620, 'Empty\n5->3')],
            'Plane 4':  [(30, 57, 'Empty\n2->3'), (665, 683, 'Empty\n3->1'), (713, 743, 'Group 61 (q=4)\n1->4')],
            'Plane 6':  [(85, 104, 'Group 7 (q=2)\n3->10')],
            'Plane 7':  [(128, 158, 'Group 21 (q=4)\n4->1'), (194, 224, 'Group 36 (q=3)\n1->4')],
            'Plane 10': [(147, 166, 'Empty\n5->4'), (420, 450, 'Group 40 (q=4)\n4->1'), (585, 603, 'Group 49 (q=2)\n1->3')],
            'Plane 12': [(187, 201, 'Group 29 (q=3)\n6->5')],
            'Plane 13': [(639, 674, 'Group 56 (q=2)\n7->5'), (704, 723, 'Empty\n5->4')],
            'Plane 15': [(186, 221, 'Group 33 (q=4)\n8->10'), (540, 575, 'Group 47 (q=2)\n10->8')],
            'Plane 17': [(252, 276, 'Empty\n9->5'), (476, 511, 'Group 42 (q=4)\n5->7')]
        }
    },
    {
        "title": "Scenario 7 (Max Profit: 13,896.14)",
        "flights": {
            'Plane 1':  [(415, 445, 'Group 31 (q=4)\n1->4'), (533, 563, 'Group 38 (q=3)\n4->1')],
            'Plane 10': [(161, 174, 'Group 24 (q=4)\n5->3')]
        }
    },
    {
        "title": "Scenario 8 (Max Profit: 49,634.44)",
        "flights": {
            'Plane 2':  [(490, 500, 'Group 40 (q=3)\n1->10')],
            'Plane 5':  [(15, 33, 'Empty\n3->1'), (100, 129, 'Group 5 (q=3)\n1->8'), (480, 509, 'Group 37 (q=4)\n8->1')],
            'Plane 6':  [(313, 326, 'Group 29 (q=4)\n3->5'), (500, 535, 'Group 42 (q=4)\n5->7'), (626, 661, 'Group 56 (q=4)\n7->5'), (676, 689, 'Empty\n5->3'), (704, 717, 'Empty\n3->5')],
            'Plane 9':  [(119, 135, 'Group 14 (q=3)\n5->2')],
            'Plane 10': [(480, 493, 'Group 38 (q=3)\n5->3'), (508, 535, 'Empty\n3->2')],
            'Plane 11': [(173, 199, 'Group 22 (q=3)\n6->3'), (483, 510, 'Group 41 (q=4)\n3->2'), (525, 545, 'Empty\n2->6'), (560, 580, 'Empty\n6->2'), (595, 611, 'Empty\n2->5')],
            'Plane 20': [(444, 463, 'Group 34 (q=3)\n10->3'), (545, 561, 'Group 46 (q=3)\n3->4'), (586, 616, 'Group 49 (q=4)\n4->1')]
        }
    },
    {
        "title": "Scenario 9 (Max Profit: 24,714.80)",
        "flights": {
            'Plane 1':  [(30, 60, 'Empty\n1->4'), (90, 108, 'Empty\n4->7'), (138, 156, 'Group 20 (q=4)\n7->1'), (186, 216, 'Empty\n4->1')],
            'Plane 2':  [(439, 469, 'Group 35 (q=4)\n1->4'), (500, 525, 'Group 41 (q=2)\n4->10'), (555, 565, 'Empty\n10->1')],
            'Plane 6':  [(54, 67, 'Group 6 (q=4)\n3->5'), (97, 110, 'Empty\n5->3'), (140, 159, 'Empty\n3->10'), (189, 208, 'Group 29 (q=4)\n10->3')],
            'Plane 7':  [(30, 45, 'Empty\n4->9'), (93, 121, 'Group 10 (q=4)\n9->3'), (734, 762, 'Empty\n3->9')],
            'Plane 9':  [(2, 21, 'Empty\n5->4'), (51, 81, 'Group 4 (q=4)\n4->1'), (111, 141, 'Group 11 (q=3)\n1->3'), (171, 184, 'Empty\n5->3')],
            'Plane 10': [(2, 32, 'Empty\n5->1'), (62, 80, 'Group 3 (q=3)\n1->3')],
            'Plane 11': [(0, 26, 'Empty\n6->3'), (56, 82, 'Empty\n3->6'), (112, 126, 'Empty\n6->5')],
            'Plane 12': [(30, 44, 'Empty\n6->5'), (74, 93, 'Empty\n5->4'), (123, 153, 'Group 18 (q=4)\n4->1'), (183, 213, 'Empty\n1->5')]
        }
    },
    {
        "title": "Scenario 10 (Max Profit: 43,672.05)",
        "flights": {
            'Plane 1':  [(175, 185, 'Group 50 (q=1)\n1->10'), (200, 210, 'Empty\n10->1'), (225, 235, 'Empty\n1->10')],
            'Plane 4':  [(704, 731, 'Group 61 (q=3)\n2->3'), (746, 769, 'Empty\n3->5')],
            'Plane 6':  [(113, 129, 'Group 8 (q=4)\n3->4'), (580, 599, 'Group 48 (q=4)\n4->1')],
            'Plane 9':  [(144, 158, 'Empty\n5->6')],
            'Plane 10': [(56, 69, 'Group 3 (q=4)\n5->3'), (181, 197, 'Group 16 (q=4)\n3->4'), (212, 228, 'Empty\n4->3'), (243, 269, 'Group 42 (q=4)\n3->10')],
            'Plane 11': [(493, 507, 'Group 40 (q=4)\n6->5'), (522, 536, 'Empty\n5->6'), (550, 564, 'Group 59 (q=3)\n6->5'), (700, 723, 'Empty\n5->6')],
            'Plane 17': [(100, 124, 'Group 5 (q=4)\n9->5'), (554, 573, 'Empty\n5->4'), (588, 603, 'Empty\n4->9')],
            'Plane 19': [(160, 187, 'Group 27 (q=4)\n10->3'), (450, 481, 'Group 38 (q=4)\n3->6'), (507, 523, 'Empty\n6->3')],
            'Plane 20': [(522, 541, 'Empty\n10->3'), (558, 575, 'Empty\n3->10'), (590, 609, 'Group 54 (q=2)\n10->4'), (626, 656, 'Group 28 (q=4)\n4->1'), (700, 725, 'Empty\n1->10')]
        }
    },
]

# ==========================================
# RENDER LOGIC
# ==========================================
num_solutions = len(all_solutions)
panel_planes = [len(s["flights"]) for s in all_solutions]
fig, axes = plt.subplots(nrows=num_solutions, ncols=1,
                         figsize=(15, 1.0 * sum(panel_planes) + 1.2 * num_solutions),
                         sharex=True,
                         gridspec_kw={'height_ratios': panel_planes})

if num_solutions == 1:
    axes = [axes]

def clean_label_text(text):
    """Transforms 'Group 34 (q=3)\n10->3' into 'Gr. 34 (3)\n10->3'"""
    if "Empty" in text:
        return text
    lines = text.split('\n')
    top_line = lines[0]
    route_line = lines[1] if len(lines) > 1 else ""
    group_num = top_line.split('Group ')[1].split(' ')[0]
    q_val = top_line.split('(q=')[1].replace(')', '')
    new_top = f"Gr. {group_num} ({q_val})"
    return f"{new_top}\n{route_line}" if route_line else new_top

def drop_decimals(title):
    """Removes the two decimals from numbers like 25,432.53 -> 25,432"""
    def repl(m):
        value = float(m.group(0).replace(',', ''))
        return f"{int(value):,}"
    return re.sub(r'\d[\d,]*\.\d+', repl, title)

for s_idx, sol in enumerate(all_solutions):
    ax = axes[s_idx]
    ax.set_axisbelow(True)
    flights = sol["flights"]

    sorted_planes = sorted(list(flights.keys()), key=lambda x: int(x.split()[-1]))

    for i, plane in enumerate(sorted_planes):
        legs = flights[plane]
        sorted_legs = sorted(legs, key=lambda x: x[0])

        placed = []  # (x_center, level) of labels already positioned

        for leg_idx, (start, end, info) in enumerate(sorted_legs):
            duration = end - start
            center_time = start + duration / 2  # label sits over the bar centre

            is_empty = "Empty" in info
            face = LIGHT_RED if is_empty else PLANE_COLORS[i % len(PLANE_COLORS)]

            ax.barh(plane, duration, left=start, color=face,
                    edgecolor=DARK_RED, lw=0.8, height=0.4,
                    hatch='////' if is_empty else None)

            short_info = clean_label_text(info)

            # Keep each label centred over its own bar. If it would overlap
            # a neighbouring label, lift it to a higher level rather than
            # shifting it sideways, and draw a connector back to the bar.
            half_w = 19.0
            level = 0
            while any(lv == level and abs(center_time - xc) < 2 * half_w
                      for xc, lv in placed):
                level += 1
            level = min(level, 2)
            placed.append((center_time, level))

            label_y = i + 0.24 + level * 0.32
            if level > 0:
                ax.plot([center_time, center_time], [i + 0.20, label_y],
                        color=DARK_RED, lw=0.5, ls=':', alpha=0.6, zorder=2)

            ax.text(
                center_time,
                label_y,
                short_info,
                ha='center',
                va='bottom',
                color=DARK_RED,
                weight='bold',
                size=7,
                zorder=4,
                bbox=dict(boxstyle="round,pad=0.15", fc="white",
                          ec=DARK_RED, lw=0.5, alpha=0.92)
            )

    ax.set_title(drop_decimals(sol["title"]), fontsize=12, pad=18, weight='bold',
                 loc='left', color=DARK_RED)
    ax.grid(axis='x', linestyle='--', alpha=0.4, color=GRID_GRAY)
    ax.set_ylim(-0.6, len(sorted_planes) + 0.4)
    ax.set_xlim(0, 800)
    for side in ('top', 'right'):
        ax.spines[side].set_visible(False)
    for side in ('left', 'bottom'):
        ax.spines[side].set_color(GRID_GRAY)

# Legend (shown once, on the first panel). Colour now encodes the plane,
# so the legend only needs to explain the deadhead hatch.
legend_handles = [
    Patch(facecolor=LIGHT_RED, edgecolor=DARK_RED, hatch='////',
          label='Deadhead (empty)')
]
axes[0].legend(handles=legend_handles, loc='upper right',
               fontsize=9, framealpha=0.95, edgecolor=GRID_GRAY)

axes[-1].set_xlabel('Time (minutes / time units)', fontsize=11,
                    labelpad=10, color=DARK_RED)
plt.tight_layout()

plt.savefig(SAVE_PATH, dpi=300, bbox_inches='tight', facecolor='white')
print(f"Themed chart saved to: {SAVE_PATH}")