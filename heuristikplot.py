import matplotlib.pyplot as plt
import numpy as np

FS       = 13
FS_TITLE = 15
FS_SUB   = 12
FS_NODE  = 12

MID_RED   = '#C95F5F'
DARK_RED  = '#641818'
GRAY      = '#888888'
BG        = 'white'

plt.rcParams.update({"font.family": "serif"})

def draw_vertiport(ax, x, y, label, color):
    ax.plot(x, y, 'o', markersize=20, color=color, zorder=5,
            markeredgecolor='white', markeredgewidth=1.5)
    ax.text(x, y, label, ha='center', va='center',
            fontsize=FS_NODE, fontweight='bold', color='white', zorder=6)

def draw_arrow(ax, p1, p2, color, lw=2.2):
    ax.annotate('', xy=p2, xytext=p1,
                arrowprops=dict(arrowstyle='->', color=color,
                                lw=lw,
                                shrinkA=10,
                                shrinkB=10,
                                connectionstyle='arc3,rad=0.0'))

def setup_ax(ax, title):
    ax.set_xlim(0, 1.0)
    ax.set_ylim(0.0, 1.0)
    ax.axis('off')
    ax.set_facecolor(BG)
    ax.set_title(title, fontsize=FS_SUB, fontweight='bold',
                 pad=6, color='#333333')

def center_arrow(fig, label):
    fig.text(0.50, 0.68, label, ha='center', va='center',
             fontsize=11, color=DARK_RED, style='italic')
    fig.text(0.50, 0.55, '→', ha='center', va='center',
             fontsize=22, color=DARK_RED, fontweight='bold')

def simple_label(ax, text, color, y=0.08):
    ax.text(0.50, y, text, ha='center', va='center',
            fontsize=FS, color=color, style='italic',
            transform=ax.transAxes)

def multicolor_label(ax, parts, y=0.08):
    char_w = 0.022
    full = ''.join(t for t, c in parts)
    x = 0.50 - len(full) * char_w / 2
    for text, color in parts:
        ax.text(x, y, text, ha='left', va='center',
                fontsize=FS, color=color, style='italic',
                transform=ax.transAxes)
        x += len(text) * char_w

###############################################################################
# Vertiport positions
###############################################################################
verts = {
    1: (0.08, 0.60),
    2: (0.33, 0.82),
    3: (0.67, 0.72),
    4: (0.92, 0.45),
}

###############################################################################
# FIGURE 1: REMOVE
###############################################################################
fig1, axes1 = plt.subplots(1, 2, figsize=(9, 3.0))

# Before
ax = axes1[0]
setup_ax(ax, 'Before')
for vid, (x, y) in verts.items():
    draw_vertiport(ax, x, y, str(vid), DARK_RED)
draw_arrow(ax, verts[1], verts[2], DARK_RED)
draw_arrow(ax, verts[2], verts[3], DARK_RED)
draw_arrow(ax, verts[3], verts[4], DARK_RED)
simple_label(ax, '1 → 2 → 3 → 4', DARK_RED)

# After
ax = axes1[1]
setup_ax(ax, 'After')
for vid, (x, y) in verts.items():
    color = GRAY if vid == 3 else DARK_RED
    draw_vertiport(ax, x, y, str(vid), color)
draw_arrow(ax, verts[1], verts[2], DARK_RED)
draw_arrow(ax, verts[2], verts[4], DARK_RED)
simple_label(ax, '1 → 2 → 4', DARK_RED)

center_arrow(fig1, 'remove V3')
fig1.text(0.50, 0.97, 'Remove Operator', ha='center', va='top',
          fontsize=FS_TITLE, fontweight='bold', color='#333333')
fig1.tight_layout(rect=[0, 0.0, 1, 0.88], w_pad=3)
plt.savefig('fig_remove.png', dpi=300, bbox_inches='tight', facecolor='white')

###############################################################################
# FIGURE 2: ADD
###############################################################################
verts_add = {
    1: (0.08, 0.60),
    2: (0.33, 0.82),
    3: (0.67, 0.72),
    4: (0.92, 0.45),
}

fig2, axes2 = plt.subplots(1, 2, figsize=(9, 3.0))

# Before
ax = axes2[0]
setup_ax(ax, 'Before')
for vid in [1, 2, 3, 4]:
    color = GRAY if vid == 2 else DARK_RED
    draw_vertiport(ax, *verts_add[vid], str(vid), color)
draw_arrow(ax, verts_add[1], verts_add[3], DARK_RED)
draw_arrow(ax, verts_add[3], verts_add[4], DARK_RED)
simple_label(ax, '1 → 3 → 4', DARK_RED)

# After
ax = axes2[1]
setup_ax(ax, 'After')
for vid, (x, y) in verts_add.items():
    color = MID_RED if vid == 2 else DARK_RED
    draw_vertiport(ax, x, y, str(vid), color)
draw_arrow(ax, verts_add[1], verts_add[2], MID_RED)
draw_arrow(ax, verts_add[2], verts_add[3], DARK_RED)
draw_arrow(ax, verts_add[3], verts_add[4], DARK_RED)
ax.text(verts_add[2][0], verts_add[2][1] + 0.11, 'new',
        ha='center', fontsize=FS, color=MID_RED, style='italic')
multicolor_label(ax, [
    ('1 → ', DARK_RED),
    ('2', MID_RED),
    (' → 3 → 4', DARK_RED),
])

center_arrow(fig2, 'insert V2')
fig2.text(0.50, 0.97, 'Add Operator', ha='center', va='top',
          fontsize=FS_TITLE, fontweight='bold', color='#333333')
fig2.tight_layout(rect=[0, 0.0, 1, 0.88], w_pad=3)
plt.savefig('fig_add.png', dpi=300, bbox_inches='tight', facecolor='white')

plt.show()
print("Figures saved.")