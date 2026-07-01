"""
Figure 2.X - Illustration of the eVTOL system and its terminology.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Ellipse

# ----- Style ----------------------------------------------------------
LIGHT = "#F9EAEA"
MID   = "#C95F5F"
DARK  = "#641818"
LABEL = "#8A3B3B"

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "Georgia"],
    "svg.fonttype": "none",
})

NODE_R   = 0.30
LEG_LW   = 2.2
SHARE_LW = 4.2


def vertiport(ax, xy, label, r=NODE_R, fs=14):
    ax.add_patch(Circle(xy, r, facecolor=LIGHT, edgecolor=DARK, lw=2.0, zorder=5))
    if label:
        ax.text(*xy, label, ha="center", va="center",
                fontsize=fs, fontweight="bold", color=DARK, zorder=6)


def leg(ax, p0, p1, lw=LEG_LW, trim0=NODE_R, trim1=NODE_R):
    p0, p1 = np.array(p0, float), np.array(p1, float)
    u = (p1 - p0) / np.hypot(*(p1 - p0))
    ax.add_patch(FancyArrowPatch(p0 + u * trim0, p1 - u * (trim1 + 0.04),
                 arrowstyle="-|>", mutation_scale=15, lw=lw, color=DARK,
                 shrinkA=0, shrinkB=0, zorder=4))


def pgroup(ax, xy, line1, line2, od, w=2.6, h=1.02):
    ax.add_patch(FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.08",
                 facecolor=LIGHT, edgecolor=MID, lw=1.4, zorder=5))
    ax.text(xy[0], xy[1] + 0.28, line1, ha="center", va="center",
            fontsize=11, fontweight="bold", color=DARK, zorder=6)
    ax.text(xy[0], xy[1] + 0.0, line2, ha="center", va="center",
            fontsize=9.5, color=DARK, zorder=6)
    ax.text(xy[0], xy[1] - 0.28, od, ha="center", va="center",
            fontsize=10, fontweight="bold", color=MID, zorder=6)


def leader(ax, p0, p1):
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=MID, lw=1.0,
            ls=(0, (3, 3)), alpha=0.85, zorder=3)


def callout(ax, xy, text, ha="left", fs=11.5):
    ax.text(*xy, text, ha=ha, va="center", fontsize=fs,
            style="italic", color=LABEL, zorder=6)


def draw_key(ax, x0, y_top):
    """Small legend explaining the generic glyphs."""
    w, h = 3.7, 3.0
    ax.add_patch(FancyBboxPatch((x0, y_top - h), w, h,
                 boxstyle="round,pad=0.04,rounding_size=0.06",
                 facecolor="white", edgecolor=MID, lw=1.2, zorder=7))
    gx, tx = x0 + 0.5, x0 + 0.95
    rows = y_top - 0.85
    step = 0.55
    ax.text(x0 + 0.25, y_top - 0.32, "Key", ha="left", va="center",
            fontsize=11.5, fontweight="bold", color=DARK, zorder=8)

    # vertiport
    ax.add_patch(Circle((gx, rows), 0.16, facecolor=LIGHT, edgecolor=DARK,
                        lw=2.0, zorder=8))
    ax.text(gx, rows, "1", ha="center", va="center",
            fontsize=9, fontweight="bold", color=DARK, zorder=9)
    ax.text(tx, rows, "vertiport", ha="left", va="center",
            fontsize=10.5, color=DARK, zorder=8)
    # flight leg
    ax.add_patch(FancyArrowPatch((gx - 0.22, rows - step), (gx + 0.22, rows - step),
                 arrowstyle="-|>", mutation_scale=12, lw=LEG_LW, color=DARK, zorder=8))
    ax.text(tx, rows - step, "flight leg", ha="left", va="center",
            fontsize=10.5, color=DARK, zorder=8)
    # shared leg
    ax.add_patch(FancyArrowPatch((gx - 0.22, rows - 2 * step), (gx + 0.22, rows - 2 * step),
                 arrowstyle="-|>", mutation_scale=12, lw=SHARE_LW, color=DARK, zorder=8))
    ax.text(tx, rows - 2 * step, "shared leg (2+ groups)", ha="left", va="center",
            fontsize=10.5, color=DARK, zorder=8)
    # passenger group
    ax.add_patch(FancyBboxPatch((gx - 0.25, rows - 3 * step - 0.16), 0.5, 0.32,
                 boxstyle="round,pad=0.01,rounding_size=0.04",
                 facecolor=LIGHT, edgecolor=MID, lw=1.3, zorder=8))
    ax.text(tx, rows - 3 * step, "passenger group", ha="left", va="center",
            fontsize=10.5, color=DARK, zorder=8)


def make_figure():
    fig, ax = plt.subplots(figsize=(11.5, 6.2))

    V1 = (2.6, 3.1)    # base
    V2 = (5.6, 4.3)
    V3 = (8.6, 3.1)
    V4 = (11.2, 4.3)   # end

    # Legs (V2->V3 shared)
    leg(ax, V1, V2)
    leg(ax, V2, V3, lw=SHARE_LW)
    leg(ax, V3, V4)

    for xy, lab in [(V1, "1"), (V2, "2"), (V3, "3"), (V4, "4")]:
        vertiport(ax, xy, lab)

    # Passenger groups below their origin nodes (wider + spaced out)
    PGA = (2.6, 1.15)
    PGB = (5.55, 1.15)
    PGC = (8.5, 1.15)
    pgroup(ax, PGA, "PG A", "2 pax, direct-only", "1 \u2192 2")
    pgroup(ax, PGB, "PG B", "3 pax, stopover allowed", "2 \u2192 4")
    pgroup(ax, PGC, "PG C", "1 pax, stopover allowed", "2 \u2192 3")
    leader(ax, (PGA[0], PGA[1] + 0.52), (V1[0], V1[1] - NODE_R))
    leader(ax, (PGB[0] - 0.3, PGB[1] + 0.52), (V2[0] - 0.15, V2[1] - NODE_R))
    leader(ax, (PGC[0] - 0.4, PGC[1] + 0.52), (V2[0] + 0.15, V2[1] - NODE_R))

    # Base vertiport (below-left of node 1, clear of the node)
    callout(ax, (0.30, 2.32), "base vertiport", fs=11)
    callout(ax, (0.30, 1.99), "(starts at 80% battery)", fs=10)
    leader(ax, (1.45, 2.25), (V1[0] - NODE_R * 0.7, V1[1] - NODE_R * 0.7))

    # End vertiport (right of node 4)
    callout(ax, (11.75, 4.46), "end-vertiport")
    callout(ax, (11.75, 4.12), "(within 60 min", fs=10)
    callout(ax, (11.75, 3.83), " drive of base)", fs=10)
    leader(ax, (V4[0] + NODE_R, V4[1]), (11.7, 4.3))

    # Stopover bracket across the top (V2 -> V3 -> V4)
    by = 5.35
    ax.plot([V2[0], V2[0], V4[0], V4[0]],
            [V2[1] + NODE_R, by, by, V4[1] + NODE_R],
            color=MID, lw=1.0, ls=(0, (3, 3)), alpha=0.8, zorder=2)
    ax.text((V2[0] + V4[0]) / 2, by + 0.12,
            "stopover trip:  2 \u2192 3 \u2192 4  (serves PG B)",
            ha="center", va="bottom", fontsize=11.5, style="italic", color=DARK)

    draw_key(ax, 0.3, 7.4)

    ax.set_xlim(0.0, 14.8)
    ax.set_ylim(0.6, 7.7)
    ax.set_aspect("equal")
    ax.axis("off")
    fig.tight_layout()
    return fig


if __name__ == "__main__":
    fig = make_figure()
    fig.savefig("system_figure.pdf", bbox_inches="tight", facecolor="white")
    fig.savefig("system_figure.png", dpi=300, bbox_inches="tight", facecolor="white")
    plt.show()
    print("Figure saved.")