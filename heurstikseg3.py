import matplotlib.pyplot as plt
import numpy as np

FS       = 11
FS_TITLE = 12
FS_NODE  = 10

MID_RED   = '#C95F5F'
DARK_RED  = '#641818'
DARK_BLUE = '#1a3a5c'
MID_BLUE  = '#7a9fc0'
GRAY      = '#888888'
BG        = 'white'

plt.rcParams.update({"font.family": "serif"})

angles = [90, 18, 306, 234, 162]
r  = 0.30
cx = 0.50

pentagon_top = {}
pentagon_bot = {}
for i, a in enumerate(angles, start=1):
    rad = np.deg2rad(a)
    pentagon_top[i] = (cx + r * np.cos(rad), 0.48 + r * np.sin(rad))
    pentagon_bot[i] = (cx + r * np.cos(rad), 0.55 + r * np.sin(rad))

def draw_vertiport(ax, x, y, label, color):
    ax.plot(x, y, 'o', markersize=16, color=color, zorder=5,
            markeredgecolor='white', markeredgewidth=1.5)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=FS_NODE, fontweight='bold', color='white', zorder=6)

def draw_arrow(ax, p1, p2, color, lw=2.0):
    ax.annotate('', xy=p2, xytext=p1,
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw,
                                shrinkA=9,
                                shrinkB=9,
                                connectionstyle='arc3,rad=0.0'))

def setup_ax(ax, title, color='#333333'):
    ax.set_xlim(0.0, 1.0)
    ax.set_ylim(0.05, 0.95)
    ax.set_aspect('equal')
    ax.axis('off')
    ax.set_facecolor(BG)
    ax.set_title(title, fontsize=FS_TITLE, fontweight='bold',
                 pad=4, color=color)

def draw_route_label(ax, segments, y=0.10):
    char_w = 0.022
    full = ''.join(t for t, c in segments)
    n = len(full)
    x = 0.50 - n * char_w / 2
    for text, color in segments:
        ax.text(x, y, text, ha='left', va='center',
                fontsize=FS, color=color, style='italic',
                transform=ax.transAxes)
        x += len(text) * char_w

def draw_pentagon_route(ax, nodes, route, normal_color,
                        highlight_legs, highlight_color):
    for vid in nodes:
        draw_vertiport(ax, *nodes[vid], str(vid), normal_color)
    for k in range(len(route) - 1):
        i, j = route[k], route[k+1]
        if (i, j) in highlight_legs:
            c, lw = highlight_color, 2.6
        else:
            c, lw = normal_color, 1.8
        draw_arrow(ax, nodes[i], nodes[j], c, lw=lw)

###############################################################################
# FIGURE: tall and narrow
###############################################################################
fig3, axes3 = plt.subplots(2, 2, figsize=(7.5, 9.0))

# Top left: eVTOL A before
ax = axes3[0, 0]
setup_ax(ax, 'eVTOL A — Before', color=DARK_RED)
draw_pentagon_route(ax, pentagon_top,
    route=[1,2,3,4,5], normal_color=DARK_RED,
    highlight_legs={(2,3),(3,4)}, highlight_color=MID_RED)
draw_route_label(ax, [
    ('1 → ', DARK_RED),
    ('2 → 3 → 4', MID_RED),
    (' → 5', DARK_RED)])

# Top right: eVTOL A after
ax = axes3[0, 1]
setup_ax(ax, 'eVTOL A — After', color=DARK_RED)
draw_pentagon_route(ax, pentagon_top,
    route=[1,4,3,2,5], normal_color=DARK_RED,
    highlight_legs={(4,3),(3,2)}, highlight_color=MID_BLUE)
draw_route_label(ax, [
    ('1 → ', DARK_RED),
    ('4 → 3 → 2', MID_BLUE),
    (' → 5', DARK_RED)])

# Bottom left: eVTOL B before
ax = axes3[1, 0]
setup_ax(ax, 'eVTOL B — Before', color=DARK_BLUE)
draw_pentagon_route(ax, pentagon_bot,
    route=[5,4,3,2,1], normal_color=DARK_BLUE,
    highlight_legs={(4,3),(3,2)}, highlight_color=MID_BLUE)
draw_route_label(ax, [
    ('5 → ', DARK_BLUE),
    ('4 → 3 → 2', MID_BLUE),
    (' → 1', DARK_BLUE)])

# Bottom right: eVTOL B after
ax = axes3[1, 1]
setup_ax(ax, 'eVTOL B — After', color=DARK_BLUE)
draw_pentagon_route(ax, pentagon_bot,
    route=[5,2,3,4,1], normal_color=DARK_BLUE,
    highlight_legs={(2,3),(3,4)}, highlight_color=MID_RED)
draw_route_label(ax, [
    ('5 → ', DARK_BLUE),
    ('2 → 3 → 4', MID_RED),
    (' → 1', DARK_BLUE)])

# Centre arrows
fig3.text(0.50, 0.735, '→', ha='center', va='center',
          fontsize=24, color=DARK_RED, fontweight='bold')
fig3.text(0.50, 0.265, '→', ha='center', va='center',
          fontsize=24, color=DARK_BLUE, fontweight='bold')

fig3.text(0.50, 0.50, 'exchange\nsegments', ha='center', va='center',
          fontsize=10, color=GRAY, style='italic',
          bbox=dict(facecolor='white', edgecolor=GRAY,
                    boxstyle='round,pad=0.3', linewidth=1.0))

fig3.text(0.50, 0.975, 'Segment Exchange Operator',
          ha='center', va='top',
          fontsize=FS_TITLE, fontweight='bold', color='#333333')

fig3.tight_layout(rect=[0, 0.0, 1, 0.955], h_pad=0.3, w_pad=0.5)
plt.savefig('fig_segment_exchange.png', dpi=300,
            bbox_inches='tight', facecolor='white')
plt.show()
print("Saved.")