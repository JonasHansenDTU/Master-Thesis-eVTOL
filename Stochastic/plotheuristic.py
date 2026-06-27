import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 13,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
})

LIGHT_RED = '#F9EAEA'
MID_RED   = '#C95F5F'
DARK_RED  = '#641818'
GRAY      = '#888888'
LGRAY     = '#cccccc'
BG        = 'white'

verts = {
    1: (0.15, 0.50),
    2: (0.50, 0.80),
    3: (0.85, 0.50),
    4: (0.50, 0.20),
}

def draw_vertiport(ax, x, y, label, color):
    ax.plot(x, y, 'o', markersize=22, color=color, zorder=5,
            markeredgecolor='white', markeredgewidth=1.5)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=12, fontweight='bold', color='white', zorder=6)

def draw_arrow(ax, p1, p2, color, lw=2.5):
    ax.annotate('', xy=p2, xytext=p1,
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw, connectionstyle='arc3,rad=0.0'))

def draw_panel(ax, title, show_stopover=False):
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_facecolor(BG)
    ax.set_title(title, fontsize=13, fontweight='bold', pad=12, color='#333333')

    if not show_stopover:
        colors = {1: DARK_RED, 2: GRAY, 3: DARK_RED, 4: GRAY}
        for vid, (x, y) in verts.items():
            draw_vertiport(ax, x, y, str(vid), colors[vid])

        draw_arrow(ax, verts[1], verts[3], color=DARK_RED)

        ax.text(verts[1][0] - 0.10, verts[1][1] + 0.09,
                'PG 1', fontsize=10, color=DARK_RED, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=LIGHT_RED,
                          edgecolor=DARK_RED, linewidth=1.2))

        ax.text(0.50, 0.58, 'Direct route', ha='center',
                fontsize=10, color=DARK_RED, style='italic')

    else:
        colors = {1: DARK_RED, 2: MID_RED, 3: DARK_RED, 4: GRAY}
        for vid, (x, y) in verts.items():
            draw_vertiport(ax, x, y, str(vid), colors[vid])

        draw_arrow(ax, verts[1], verts[2], color=MID_RED)
        draw_arrow(ax, verts[2], verts[3], color=MID_RED)

        ax.text(verts[1][0] - 0.10, verts[1][1] + 0.09,
                'PG 1', fontsize=10, color=DARK_RED, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=LIGHT_RED,
                          edgecolor=DARK_RED, linewidth=1.2))

        ax.text(verts[2][0] + 0.11, verts[2][1] + 0.07,
                'PG 2', fontsize=10, color=MID_RED, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor=LIGHT_RED,
                          edgecolor=MID_RED, linewidth=1.2))

        ax.text(verts[2][0], verts[2][1] - 0.12,
                'Stopover', ha='center', fontsize=9,
                color=MID_RED, style='italic')

        ax.text(0.50, 0.05, 'Shared route with stopover', ha='center',
                fontsize=10, color=MID_RED, style='italic')

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

draw_panel(axes[0], 'Option A: Direct Route', show_stopover=False)
draw_panel(axes[1], 'Option B: Shared Route with Stopover', show_stopover=True)

legend_elements = [
    mpatches.Patch(facecolor=DARK_RED, label='Direct route (PG 1)'),
    mpatches.Patch(facecolor=MID_RED,  label='Shared route (PG 1 + PG 2)'),
    mpatches.Patch(facecolor=GRAY,     label='Inactive vertiport'),
]
fig.legend(handles=legend_elements, loc='lower center', ncol=3,
           fontsize=10, frameon=True, facecolor='white',
           edgecolor=LGRAY, bbox_to_anchor=(0.5, -0.02))

plt.suptitle('Candidate Route Options for Passenger Group 1',
             fontsize=14, fontweight='bold', color='#333333', y=1.01)

fig.tight_layout()
plt.savefig('heuristic_figure.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.show()