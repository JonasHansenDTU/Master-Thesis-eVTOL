"""
Figure 2.X — Illustration of the eVTOL system and its terminology.

One eVTOL route over four vertiports, annotating the core terms used
throughout the thesis: vertiport (node), flight leg (edge), passenger
group (box), base vertiport, end-vertiport, shared leg and stopover trip.

Style matches the thesis palette and the Chapter 6 neighbourhood diagrams.
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Circle, Ellipse
from matplotlib import font_manager

# ----------------------------------------------------------------------
# Style
# ----------------------------------------------------------------------
LIGHT = "#F9EAEA"   # light fill
MID   = "#C95F5F"   # mid / accents / leaders
DARK  = "#641818"   # nodes, legs, text
LABEL = "#8A3B3B"   # italic callout text

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "Georgia"],
    "svg.fonttype": "none",
})

NODE_R   = 0.34     # vertiport radius (data units)
LEG_LW   = 2.2      # normal flight leg line width
SHARE_LW = 4.0      # shared leg line width


# ----------------------------------------------------------------------
# Helpers
# ----------------------------------------------------------------------
def vertiport(ax, xy, label):
    """Draw a numbered vertiport node."""
    ax.add_patch(Circle(xy, NODE_R, facecolor=LIGHT, edgecolor=DARK,
                        lw=2.0, zorder=5))
    ax.text(*xy, label, ha="center", va="center",
            fontsize=15, fontweight="bold", color=DARK, zorder=6)


def leg(ax, p0, p1, lw=LEG_LW):
    """Draw a flight leg as an arrow, trimmed to the node edges."""
    import numpy as np
    p0, p1 = np.array(p0, float), np.array(p1, float)
    d = p1 - p0
    u = d / np.hypot(*d)
    start = p0 + u * NODE_R
    end   = p1 - u * (NODE_R + 0.04)
    ax.add_patch(FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=16,
        lw=lw, color=DARK, shrinkA=0, shrinkB=0, zorder=4))


def pgroup(ax, xy, text, w=2.05, h=0.62):
    """Draw a passenger-group box centred at xy."""
    box = FancyBboxPatch(
        (xy[0] - w / 2, xy[1] - h / 2), w, h,
        boxstyle="round,pad=0.02,rounding_size=0.10",
        facecolor=LIGHT, edgecolor=MID, lw=1.4, zorder=5)
    ax.add_patch(box)
    ax.text(*xy, text, ha="center", va="center",
            fontsize=11.5, color=DARK, zorder=6)


def leader(ax, p0, p1):
    """Thin dashed leader line."""
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]],
            color=MID, lw=1.0, ls=(0, (3, 3)), alpha=0.8, zorder=3)


def callout(ax, xy, text, ha="left"):
    ax.text(*xy, text, ha=ha, va="center", fontsize=11.5,
            style="italic", color=LABEL, zorder=6)


# ----------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------
def make_figure():
    fig, ax = plt.subplots(figsize=(9.2, 5.6))

    # Node positions
    V1 = (1.10, 3.30)   # base
    V2 = (3.60, 4.85)
    V3 = (6.70, 3.30)
    V4 = (8.80, 4.85)   # end-vertiport

    # Drive-radius hint around the end-vertiport
    ax.add_patch(Ellipse(V4, 1.25, 1.10, facecolor="none",
                        edgecolor=MID, lw=1.0, ls=(0, (4, 4)),
                        alpha=0.7, zorder=2))

    # Flight legs (V2->V3 is shared, hence thicker)
    leg(ax, V1, V2)
    leg(ax, V2, V3, lw=SHARE_LW)
    leg(ax, V3, V4)

    # Vertiports
    vertiport(ax, V1, "1")
    vertiport(ax, V2, "2")
    vertiport(ax, V3, "3")
    vertiport(ax, V4, "4")

    # Passenger groups
    PGA = (1.10, 1.35)
    PGB = (3.60, 6.55)
    PGC = (3.60, 1.35)
    pgroup(ax, PGA, "PG A · 2 pax · direct-only")
    pgroup(ax, PGB, "PG B · 3 pax · stopover")
    pgroup(ax, PGC, "PG C · 1 pax", w=1.55)
    leader(ax, (PGA[0], PGA[1] + 0.31), (V1[0], V1[1] - NODE_R))
    leader(ax, (PGB[0], PGB[1] - 0.31), (V2[0], V2[1] + NODE_R))
    leader(ax, (PGC[0], PGC[1] + 0.31), (V2[0], V2[1] - NODE_R))

    # Callouts
    callout(ax, (0.10, 2.55), "base vertiport")
    callout(ax, (0.10, 2.25), "(starts at 80% battery)")
    leader(ax, (1.55, 2.55), (V1[0] + 0.22, V1[1] - NODE_R + 0.05))

    callout(ax, (1.95, 4.55), "flight leg")
    leader(ax, (2.35, 4.45), (2.35, 4.06))

    callout(ax, (4.95, 2.05), "shared leg (PG B + PG C)", ha="center")
    leader(ax, (4.95, 2.30), (5.15, 3.70))

    callout(ax, (6.05, 2.30), "vertiport = node", ha="center")
    leader(ax, (6.45, 2.55), (6.55, 2.98))

    callout(ax, (9.55, 4.85), "end-vertiport", ha="left")
    callout(ax, (9.55, 4.55), "(within 60 min", ha="left")
    callout(ax, (9.55, 4.28), " drive of base)", ha="left")
    leader(ax, (9.50, 4.65), (V4[0] + NODE_R, V4[1]))

    # Stopover bracket spanning V2 -> V3 -> V4
    bx0, bx1, by = V2[0], V4[0], 6.10
    ax.plot([bx0, bx0, bx1, bx1], [V2[1] + NODE_R + 0.05, by, by, V4[1] + NODE_R + 0.05],
            color=MID, lw=1.0, ls=(0, (3, 3)), alpha=0.7, zorder=2)
    ax.text((V3[0] + bx1) / 2, by + 0.20,
            "stopover trip: 2 → 3 → 4 (serves PG B)",
            ha="center", va="bottom", fontsize=12, style="italic", color=DARK)

    # Frame
    ax.set_xlim(-0.2, 11.6)
    ax.set_ylim(0.6, 7.0)
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