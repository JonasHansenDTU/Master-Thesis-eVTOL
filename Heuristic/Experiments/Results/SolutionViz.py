import os
import matplotlib.pyplot as plt

# Get the directory where this script is saved
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(SCRIPT_DIR, 'evtol_schedule_english.png')

# ==========================================
# MASTER DATA CONFIGURATION
# ==========================================
all_solutions = [
    # Solution 1
    {
        "title": "Base case (Max Profit: 25,432.53)",
        "flights": {
            'Plane 3':  [(153, 180, 'Group 21 (q=4)\n2->3'), (611, 638, 'Group 52 (q=4)\n3->2')],
            'Plane 5':  [(106, 119, 'Group 7 (q=1)\n3->5'),  (149, 162, 'Group 18 (q=3)\n5->3')],
            'Plane 10': [(216, 230, 'Group 27 (q=3)\n5->6'), (260, 274, 'Empty\n6->5')],
            'Plane 20': [(92, 117, 'Group 2 (q=4)\n10->4'),  (722, 752, 'Group 60 (q=4)\n4->1')]
        }
    },
    # Solution 2
    {
        "title": "Scenario 1 (Max Profit: 31,820.06)",
        "flights": {
            'Plane 1':  [(105, 115, 'Group 7 (q=3)\n1->10'), (622, 641, 'Empty\n10->3'), (656, 674, 'Empty\n3->1')],
            'Plane 5':  [(395, 411, 'Group 26 (q=3)\n3->4'), (477, 507, 'Group 29 (q=4)\n4->1'), (522, 540, 'Empty\n1->3')],
            'Plane 19': [(15, 25, 'Group 46 (q=3)\n10->3'), (79, 109, 'Group 5 (q=3)\n1->4'), (494, 524, 'Group 30 (q=4)\n4->1'), (560, 590, 'Group 39 (q=4)\n1->4'), (689, 719, 'Group 54 (q=4)\n4->1')]
        }
    },
    # Solution 3
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
    # Solution 4
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
    # Scenario 4
    {
        "title": "Scenario 4 (Max Profit: 18,041.75)",
        "flights": {
            'Plane 7':  [(42, 72, 'Group 3 (q=4)\n4->1'), (102, 132, 'Group 14 (q=4)\n1->4')],
            'Plane 8':  [(127, 157, 'Group 21 (q=3)\n4->1'), (358, 388, 'Group 30 (q=2)\n1->4')]
        }
    },
    # Solution 5
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
    # Solution 6
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
    # Scenario 7
    {
        "title": "Scenario 7 (Max Profit: 13,896.14)",
        "flights": {
            'Plane 1':  [(415, 445, 'Group 31 (q=4)\n1->4'), (533, 563, 'Group 38 (q=3)\n4->1')],
            'Plane 10': [(161, 174, 'Group 24 (q=4)\n5->3')]
        }
    },
    # Solution 7
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
    # Scenario 9
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
    # Scenario 10
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
fig, axes = plt.subplots(nrows=num_solutions, ncols=1, figsize=(15, 4.5 * num_solutions), sharex=True)

if num_solutions == 1:
    axes = [axes]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']

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

for s_idx, sol in enumerate(all_solutions):
    ax = axes[s_idx]
    flights = sol["flights"]
    
    sorted_planes = sorted(list(flights.keys()), key=lambda x: int(x.split()[-1]))
    
    for i, plane in enumerate(sorted_planes):
        legs = flights[plane]
        sorted_legs = sorted(legs, key=lambda x: x[0])
        
        last_text_end = -999.0
        
        for leg_idx, (start, end, info) in enumerate(sorted_legs):
            duration = end - start
            center_time = start + duration / 2
            
            ax.barh(plane, duration, left=start, color=colors[i % len(colors)], edgecolor='black', alpha=0.85, height=0.4)
            
            short_info = clean_label_text(info)
            
            # Anti-collision layout rules
            estimated_box_half_width = 19.0
            if center_time - estimated_box_half_width < last_text_end:
                center_time = last_text_end + estimated_box_half_width + 4.0
                
            last_text_end = center_time + estimated_box_half_width
            
            # Text boxes locked strictly above the bars
            ax.text(
                center_time, 
                i + 0.28,           
                short_info, 
                ha='center', 
                va='bottom', 
                color='black',      
                weight='bold',
                size=7,
                bbox=dict(boxstyle="round,pad=0.15", fc="white", ec="gray", lw=0.5, alpha=0.92)
            )
            
    ax.set_title(sol["title"], fontsize=12, pad=18, weight='bold', loc='left')
    ax.grid(axis='x', linestyle='--', alpha=0.4)
    ax.set_ylim(-0.6, len(sorted_planes) - 0.2)
    ax.set_xlim(0, 800)

axes[-1].set_xlabel('Time (minutes / time units)', fontsize=11, labelpad=10)
plt.tight_layout()

plt.savefig(SAVE_PATH, dpi=300, bbox_inches='tight')
print(f"Clean chart updated and saved directly to: {SAVE_PATH}")