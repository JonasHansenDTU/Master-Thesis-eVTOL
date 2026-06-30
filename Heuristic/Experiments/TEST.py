import os
import matplotlib.pyplot as plt

# Get the directory where this script is saved
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
SAVE_PATH = os.path.join(SCRIPT_DIR, 'single_evtol_schedule.png')

# ==========================================
# SOLUTION DATA CONFIGURATION
# ==========================================
solution = {
    "title": "Optimal Solution (Max Profit: 45,532.95)",
    "flights": {
        'Plane 1':  [(192, 222, 'Group 15 (q=3)\n1->4'), (497, 527, 'Group 24 (q=3)\n4->1')],
        'Plane 4':  [(206, 222, 'Empty\n2->5'), (669, 683, 'Group 50 (q=2)\n5->6')],
        'Plane 5':  [(77, 95, 'Empty\n3->1'), (303, 333, 'Group 18 (q=4)\n1->4'), (497, 527, 'Empty\n4->1'), (646, 676, 'Group 47 (q=4)\n1->4'), (716, 746, 'Group 55 (q=4)\n4->1')],
        'Plane 12': [(76, 109, 'Empty\n6->4'), (160, 176, 'Group 7 (q=4)\n4->3'), (191, 209, 'Group 7 (q=4)\n3->1'), (550, 580, 'Group 32 (q=4)\n1->4'), (595, 614, 'Empty\n4->5')],
        'Plane 13': [(171, 206, 'Group 12 (q=4)\n7->5'), (511, 525, 'Group 27 (q=2)\n5->6'), (555, 588, 'Empty\n6->4')],
        'Plane 19': [(15, 40, 'Empty\n10->4'), (71, 101, 'Group 2 (q=3)\n4->1')]
    }
}

# ==========================================
# RENDER LOGIC
# ==========================================
fig, ax = plt.subplots(figsize=(15, 6.0))

colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']

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

flights = solution["flights"]
sorted_planes = sorted(list(flights.keys()), key=lambda x: int(x.split()[-1]))

for i, plane in enumerate(sorted_planes):
    legs = flights[plane]
    sorted_legs = sorted(legs, key=lambda x: x[0])
    
    last_text_end = -999.0
    
    for leg_idx, (start, end, info) in enumerate(sorted_legs):
        duration = end - start
        center_time = start + duration / 2
        
        # Plot aircraft active timeline bar
        ax.barh(plane, duration, left=start, color=colors[i % len(colors)], edgecolor='black', alpha=0.85, height=0.4)
        
        short_info = clean_label_text(info)
        
        # Text alignment collision guide
        estimated_box_half_width = 18.0
        if center_time - estimated_box_half_width < last_text_end:
            center_time = last_text_end + estimated_box_half_width + 4.0
            
        last_text_end = center_time + estimated_box_half_width
        
        # Display readable data labels above the schedule bars
        ax.text(
            center_time, 
            i + 0.26,           
            short_info, 
            ha='center', 
            va='bottom', 
            color='black',      
            weight='bold',
            size=7.5,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="gray", lw=0.5, alpha=0.95)
        )
        
ax.set_title(solution["title"], fontsize=14, pad=20, weight='bold', loc='left')
ax.grid(axis='x', linestyle='--', alpha=0.4)
ax.set_ylim(-0.6, len(sorted_planes) - 0.2)
ax.set_xlim(0, 800)
ax.set_xlabel('Time (minutes / time units)', fontsize=11, labelpad=12)

plt.tight_layout()

# Save the visualization directly
plt.savefig(SAVE_PATH, dpi=300, bbox_inches='tight')
print(f"Clean timeline schedule chart saved to: {SAVE_PATH}")