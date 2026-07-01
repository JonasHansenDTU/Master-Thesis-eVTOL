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
        "title": "Base case (Max Profit: 32,507.99)",
        "flights": {
            'Plane 2':  [(489, 519, 'Group 32 (q=3)\n1->4'), (555, 585, 'Group 43 (q=4)\n4->1'), (615, 645, 'Group 49 (q=3)\n1->4'), (675, 705, 'Group 56 (q=2)\n4->1')],
            'Plane 6':  [(0, 18, 'Empty\n3->1'), (198, 228, 'Group 20 (q=4)\n1->4'), (473, 489, 'Group 29 (q=4)\n4->3')],
            'Plane 7':  [(173, 203, 'Group 14 (q=4)\n4->1'), (388, 418, 'Group 23 (q=4)\n1->4')]
        }
    },
    # {
    #     "title": "Scenario 1 (Max Profit: 31,649.83)",
    #     "flights": {
    #         'Plane 1':  [(151, 180, 'Group 15 (q=4)\n1->8'), (590, 619, 'Group 47 (q=4)\n8->1')],
    #         'Plane 6':  [(22, 38, 'Empty\n3->4'), (200, 230, 'Group 25 (q=4)\n4->1'), (245, 255, 'Empty\n1->10'), (501, 520, 'Group 37 (q=4)\n10->3')],
    #         'Plane 8':  [(111, 127, 'Group 11 (q=4)\n4->3'), (142, 160, 'Group 11 (q=4)\n3->1'), (598, 628, 'Group 46 (q=4)\n1->4')],
    #         'Plane 10': [(75, 88, 'Empty\n5->3'), (570, 596, 'Group 44 (q=3)\n3->6'), (613, 627, 'Empty\n6->5')],
    #         'Plane 16': [(81, 116, 'Group 7 (q=4)\n8->10'), (422, 457, 'Group 30 (q=3)\n10->8')]
    #     }
    # },
    # {
    #     "title": "Scenario 2 (Max Profit: 61,597.75)",
    #     "flights": {
    #         'Plane 1':  [(130, 160, 'Group 8 (q=4)\n1->4'), (588, 618, 'Group 41 (q=2)\n4->1')],
    #         'Plane 5':  [(36, 54, 'Empty\n3->1'), (133, 163, 'Group 10 (q=3)\n1->4'), (193, 216, 'Empty\n4->2'), (597, 624, 'Group 45 (q=3)\n2->3')],
    #         'Plane 7':  [(146, 176, 'Group 18 (q=4)\n4->1'), (543, 573, 'Group 34 (q=4)\n1->4')],
    #         'Plane 8':  [(0, 15, 'Empty\n4->9'), (613, 653, 'Group 47 (q=2)\n9->10'), (706, 731, 'Group 56 (q=3)\n10->4')],
    #         'Plane 9':  [(92, 123, 'Group 5 (q=4)\n5->10'), (564, 574, 'Group 38 (q=4)\n10->1'), (738, 768, 'Empty\n1->5')],
    #         'Plane 11': [(2, 28, 'Empty\n6->3'), (422, 448, 'Group 29 (q=4)\n3->6')],
    #         'Plane 13': [(125, 141, 'Empty\n7->9'), (171, 195, 'Group 19 (q=4)\n9->5'), (544, 574, 'Empty\n5->1'), (604, 634, 'Group 42 (q=3)\n1->4'), (664, 682, 'Group 42 (q=3)\n4->7')]
    #     }
    # },
    # {
    #     "title": "Scenario 3 (Max Profit: 33,321.68)",
    #     "flights": {
    #         'Plane 1':  [(167, 197, 'Group 21 (q=3)\n1->4'), (270, 285, 'Empty\n4->9'), (507, 547, 'Group 33 (q=4)\n9->10')],
    #         'Plane 2':  [(159, 189, 'Group 19 (q=4)\n1->4'), (610, 640, 'Group 46 (q=4)\n4->1')],
    #         'Plane 5':  [(531, 549, 'Group 35 (q=4)\n3->1'), (580, 609, 'Group 40 (q=4)\n1->8'), (738, 767, 'Empty\n8->1')],
    #         'Plane 20': [(113, 148, 'Group 12 (q=4)\n10->8'), (178, 207, 'Group 20 (q=4)\n8->1')]
    #     }
    # },
    # {
    #     "title": "Scenario 4 (Max Profit: 32,327.73)",
    #     "flights": {
    #         'Plane 1':  [(60, 90, 'Group 3 (q=4)\n1->4'), (657, 687, 'Group 50 (q=4)\n4->1')],
    #         'Plane 2':  [(250, 280, 'Group 23 (q=4)\n1->4'), (545, 575, 'Group 39 (q=4)\n4->1')],
    #         'Plane 4':  [(40, 60, 'Empty\n2->6'), (124, 150, 'Group 12 (q=3)\n6->3'), (565, 592, 'Group 40 (q=4)\n3->2')],
    #         'Plane 15': [(540, 569, 'Group 37 (q=4)\n8->1'), (618, 647, 'Empty\n1->8')],
    #         'Plane 17': [(68, 96, 'Group 5 (q=4)\n9->3'), (231, 259, 'Empty\n3->9')]
    #     }
    # },
    # {
    #     "title": "Scenario 5 (Max Profit: 68,477.61)",
    #     "flights": {
    #         'Plane 1':  [(77, 106, 'Group 4 (q=4)\n1->8'), (476, 505, 'Group 35 (q=4)\n8->1')],
    #         'Plane 3':  [(47, 74, 'Empty\n2->3'), (463, 482, 'Group 33 (q=4)\n3->10'), (673, 713, 'Group 54 (q=4)\n10->9'), (728, 744, 'Empty\n9->2')],
    #         'Plane 6':  [(42, 60, 'Empty\n3->1'), (75, 105, 'Group 2 (q=4)\n1->4'), (120, 136, 'Empty\n4->3')],
    #         'Plane 7':  [(99, 129, 'Group 8 (q=4)\n4->1'), (162, 192, 'Group 19 (q=3)\n1->4')],
    #         'Plane 17': [(72, 87, 'Empty\n9->4'), (102, 132, 'Group 7 (q=3)\n4->1'), (147, 157, 'Empty\n1->10'), (348, 388, 'Empty\n10->9')],
    #         'Plane 19': [(160, 179, 'Group 17 (q=4)\n10->3'), (194, 210, 'Empty\n3->4'), (225, 241, 'Group 24 (q=4)\n4->3'), (256, 274, 'Group 24 (q=4)\n3->1'), (584, 614, 'Group 42 (q=4)\n1->4'), (672, 702, 'Group 52 (q=4)\n4->1')]
    #     }
    # },
    # {
    #     "title": "Scenario 6 (Max Profit: 64,528.20)",
    #     "flights": {
    #         'Plane 2':  [(166, 196, 'Group 21 (q=4)\n1->4'), (521, 551, 'Group 37 (q=3)\n4->1')],
    #         'Plane 5':  [(88, 104, 'Empty\n3->4'), (150, 180, 'Group 15 (q=3)\n4->1')],
    #         'Plane 7':  [(41, 71, 'Group 3 (q=4)\n4->1'), (349, 367, 'Empty\n1->3'), (432, 460, 'Group 31 (q=4)\n3->9')],
    #         'Plane 8':  [(30, 55, 'Empty\n4->10'), (164, 189, 'Group 19 (q=3)\n10->4'), (219, 234, 'Group 19 (q=3)\n4->9')],
    #         'Plane 10': [(5, 21, 'Group 1 (q=3)\n5->2')],
    #         'Plane 14': [(271, 289, 'Empty\n7->4'), (576, 606, 'Group 42 (q=3)\n4->1'), (748, 778, 'Group 58 (q=3)\n1->4')],
    #         'Plane 15': [(574, 609, 'Empty\n8->10'), (647, 682, 'Group 48 (q=2)\n10->8')],
    #         'Plane 16': [(643, 672, 'Group 47 (q=4)\n8->1'), (725, 754, 'Empty\n1->8')],
    #         'Plane 19': [(182, 192, 'Empty\n10->1'), (245, 275, 'Group 25 (q=4)\n1->5'), (305, 336, 'Empty\n5->10'), (509, 528, 'Group 35 (q=4)\n10->3')]
    #     }
    # },
    # {
    #     "title": "Scenario 7 (Max Profit: 47,458.39)",
    #     "flights": {
    #         'Plane 1':  [(586, 616, 'Group 43 (q=4)\n1->4'), (646, 676, 'Group 48 (q=3)\n4->1')],
    #         'Plane 4':  [(116, 139, 'Empty\n2->4'), (598, 628, 'Group 45 (q=4)\n4->1'), (658, 688, 'Group 47 (q=4)\n1->4')],
    #         'Plane 6':  [(121, 137, 'Empty\n3->4'), (508, 538, 'Group 37 (q=4)\n4->1')],
    #         'Plane 7':  [(40, 70, 'Group 1 (q=4)\n4->1'), (102, 132, 'Group 8 (q=4)\n1->4')],
    #         'Plane 15': [(81, 110, 'Group 7 (q=4)\n8->1'), (739, 768, 'Group 56 (q=4)\n1->8')],
    #         'Plane 16': [(575, 610, 'Group 42 (q=4)\n8->10'), (640, 675, 'Empty\n10->8')]
    #     }
    # },
    # {
    #     "title": "Scenario 8 (Max Profit: 82,354.61)",
    #     "flights": {
    #         'Plane 1':  [(61, 91, 'Group 3 (q=4)\n1->4'), (310, 325, 'Empty\n4->9'), (657, 697, 'Group 56 (q=4)\n9->10')],
    #         'Plane 2':  [(488, 518, 'Group 43 (q=3)\n1->4'), (597, 627, 'Group 50 (q=3)\n4->1'), (648, 678, 'Group 54 (q=3)\n1->4'), (693, 709, 'Empty\n4->3')],
    #         'Plane 6':  [(75, 91, 'Empty\n3->4'), (610, 626, 'Group 52 (q=4)\n4->3')],
    #         'Plane 8':  [(53, 83, 'Group 1 (q=4)\n4->1'), (106, 136, 'Group 14 (q=1)\n1->4')],
    #         'Plane 10': [(32, 63, 'Empty\n5->10'), (119, 150, 'Group 13 (q=4)\n10->5')],
    #         'Plane 11': [(0, 20, 'Empty\n6->2'), (116, 132, 'Group 17 (q=1)\n2->5')],
    #         'Plane 12': [(0, 14, 'Empty\n6->5'), (135, 151, 'Group 24 (q=3)\n5->2')],
    #         'Plane 13': [(120, 155, 'Group 19 (q=4)\n7->5'), (354, 373, 'Empty\n5->4'), (472, 502, 'Group 41 (q=3)\n4->1'), (598, 628, 'Group 51 (q=4)\n1->4')],
    #         'Plane 17': [(49, 64, 'Empty\n9->4'), (79, 109, 'Group 2 (q=4)\n4->1'), (133, 163, 'Group 23 (q=4)\n1->4')]
    #     }
    # },
    # {
    #     "title": "Scenario 9 (Max Profit: 66,314.96)",
    #     "flights": {
    #         'Plane 2':  [(91, 121, 'Group 12 (q=4)\n1->4'), (519, 549, 'Group 36 (q=4)\n4->1')],
    #         'Plane 3':  [(38, 61, 'Empty\n2->4'), (112, 142, 'Group 17 (q=3)\n4->1'), (640, 670, 'Group 51 (q=4)\n1->4')],
    #         'Plane 7':  [(85, 118, 'Empty\n4->6'), (151, 165, 'Group 20 (q=4)\n6->5'), (605, 636, 'Group 48 (q=4)\n5->10'), (667, 686, 'Empty\n10->3'), (716, 729, 'Empty\n3->5'), (759, 775, 'Empty\n5->2')],
    #         'Plane 10': [(30, 65, 'Group 1 (q=4)\n5->7'), (523, 541, 'Empty\n7->4'), (571, 601, 'Group 39 (q=4)\n4->1'), (631, 661, 'Group 50 (q=3)\n1->4'), (691, 710, 'Empty\n4->5')],
    #         'Plane 15': [(563, 598, 'Group 42 (q=4)\n8->10'), (662, 697, 'Empty\n10->8')],
    #         'Plane 17': [(534, 549, 'Empty\n9->4'), (602, 632, 'Group 47 (q=3)\n4->1'), (662, 692, 'Group 55 (q=3)\n1->4')],
    #         'Plane 18': [(106, 146, 'Empty\n9->10'), (551, 582, 'Group 38 (q=3)\n10->5'), (653, 683, 'Empty\n5->1'), (725, 755, 'Group 58 (q=3)\n1->4')]
    #     }
    # },
    # {
    #     "title": "Scenario 10 (Max Profit: 62,625.95)",
    #     "flights": {
    #         'Plane 1':  [(114, 144, 'Group 11 (q=4)\n1->4'), (358, 388, 'Group 32 (q=4)\n4->1')],
    #         'Plane 3':  [(587, 614, 'Group 46 (q=4)\n2->3'), (629, 647, 'Empty\n3->1'), (740, 770, 'Group 60 (q=4)\n1->4')],
    #         'Plane 4':  [(92, 108, 'Group 6 (q=4)\n2->5')],
    #         'Plane 8':  [(445, 461, 'Group 34 (q=4)\n4->3'), (476, 494, 'Group 34 (q=4)\n3->1'), (593, 623, 'Group 48 (q=4)\n1->4'), (638, 656, 'Group 48 (q=4)\n4->7')],
    #         'Plane 11': [(87, 113, 'Group 4 (q=1)\n6->3'), (170, 183, 'Empty\n3->5'), (635, 666, 'Group 51 (q=4)\n5->10'), (723, 754, 'Empty\n10->5')],
    #         'Plane 12': [(3, 17, 'Empty\n6->5'), (112, 128, 'Group 8 (q=3)\n5->2'), (143, 170, 'Group 19 (q=3)\n2->3'), (185, 203, 'Group 19 (q=3)\n3->1'), (218, 228, 'Empty\n1->10'), (660, 691, 'Group 56 (q=4)\n10->5')],
    #         'Plane 16': [(15, 44, 'Empty\n8->1'), (621, 650, 'Group 49 (q=3)\n1->8')],
    #         'Plane 17': [(550, 578, 'Group 43 (q=4)\n9->3'), (593, 612, 'Group 43 (q=4)\n3->10'), (752, 777, 'Empty\n10->4')],
    #         'Plane 20': [(15, 34, 'Empty\n10->3'), (104, 123, 'Group 7 (q=3)\n3->10')]
    #     }
    # },
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