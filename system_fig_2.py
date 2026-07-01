"""Figure 2.X - eVTOL system and terminology (Chapter-6-consistent style)."""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle

LIGHT = "#F9EAEA"
MID   = "#C95F5F"
DARK  = "#641818"
LABEL = "#8A3B3B"
GRAY  = "#888888"

FS_NODE = 12
FS      = 12
plt.rcParams.update({"font.family": "serif",
                     "font.serif": ["DejaVu Serif", "Times New Roman", "Georgia"],
                     "svg.fonttype": "none"})

MS = 26   # node marker size (points), matches Ch.6 style


def node(ax, xy, label, color=DARK, ms=MS, fs=FS_NODE):
    ax.plot(*xy, "o", markersize=ms, color=color, zorder=5,
            markeredgecolor="white", markeredgewidth=1.6)
    if label:
        ax.text(*xy, label, ha="center", va="center", fontsize=fs,
                fontweight="bold", color="white", zorder=6)


def arrow(ax, p0, p1, color=DARK, lw=2.2):
    ax.annotate("", xy=p1, xytext=p0, zorder=4,
                arrowprops=dict(arrowstyle="-|>", color=color, lw=lw,
                                shrinkA=13, shrinkB=13,
                                mutation_scale=16))


def pgroup(ax, xy, name, desc, od, w=2.7, h=1.0):
    ax.add_patch(FancyBboxPatch((xy[0] - w / 2, xy[1] - h / 2), w, h,
                 boxstyle="round,pad=0.02,rounding_size=0.08",
                 facecolor=LIGHT, edgecolor=MID, lw=1.4, zorder=6))
    ax.text(xy[0], xy[1] + 0.27, name, ha="center", va="center",
            fontsize=11, fontweight="bold", color=DARK, zorder=7)
    ax.text(xy[0], xy[1] - 0.01, desc, ha="center", va="center",
            fontsize=9.3, color=DARK, zorder=7)
    ax.text(xy[0], xy[1] - 0.28, od, ha="center", va="center",
            fontsize=10, fontweight="bold", color=MID, zorder=7)


def leader(ax, p0, p1):
    ax.plot([p0[0], p1[0]], [p0[1], p1[1]], color=MID, lw=1.0,
            ls=(0, (3, 3)), alpha=0.85, zorder=3)


def callout(ax, xy, text, ha="center", fs=11, color=LABEL):
    ax.text(*xy, text, ha=ha, va="center", fontsize=fs, style="italic",
            color=color, zorder=8)


def draw_legend(ax, x0, y):
    """Horizontal key strip along the bottom."""
    w, h = 13.1, 0.5
    ax.add_patch(FancyBboxPatch((x0, y - h), w, h,
                 boxstyle="round,pad=0.05,rounding_size=0.06",
                 facecolor="white", edgecolor=MID, lw=1.2, zorder=9))
    cy = y - h / 2
    ax.text(x0 + 0.45, cy, "Key", ha="left", va="center",
            fontsize=11.5, fontweight="bold", color=DARK, zorder=11)

    def marker(px, color):
        ax.plot(px, cy, "o", markersize=15, color=color,
                markeredgecolor="white", markeredgewidth=1.4, zorder=11)

    def arr(px, lw):
        ax.annotate("", xy=(px + 0.32, cy), xytext=(px - 0.32, cy),
                    arrowprops=dict(arrowstyle="-|>", color=DARK, lw=lw,
                                    mutation_scale=13), zorder=11)

    def txt(px, s):
        ax.text(px, cy, s, ha="left", va="center", fontsize=10,
                color=DARK, zorder=11)

    marker(x0 + 1.55, DARK);  txt(x0 + 1.8, "vertiport")
    marker(x0 + 3.25, GRAY);  txt(x0 + 3.5, "alternative end-vertiport")
    arr(x0 + 6.60, 1.45);      txt(x0 + 7.05, "flight leg")
    arr(x0 + 8.65, 3.10);      txt(x0 + 9.10, "shared leg")
    ax.add_patch(FancyBboxPatch((x0 + 10.75, cy - 0.16), 0.5, 0.32,
                 boxstyle="round,pad=0.01,rounding_size=0.04",
                 facecolor=LIGHT, edgecolor=MID, lw=1.2, zorder=11))
    txt(x0 + 11.40, "passenger group")


def make_figure():
    fig, ax = plt.subplots(figsize=(12.6, 7.4))

    V1 = (3.5, 4.0)   # Aalborg
    V2 = (8.8, 1.9)   # Kalundborg
    V3 = (6.0, 1.0)   # Odense
    V4 = (4.2, 2.6)   # Aarhus

    # alternative end-vertiports
    G1 = (2.0, 3.0)   # sydvest for Aalborg
    G2 = (2.2, 4.9)   # nordvest for Aalborg
    R  = 2.05

    # 60-minute drive region, centred on the base
    ax.add_patch(Circle(V1, R, edgecolor=MID, lw=1.2, ls=(0, (5, 4)),
                        fill=False, zorder=2))
    ax.add_patch(Circle(V1, R, facecolor=MID, alpha=0.06, zorder=1))

    # Legs: outbound, shared (vertical), return -> symmetric about y = 3.4
    arrow(ax, V1, V2)
    arrow(ax, V2, V3, lw=4.2)
    arrow(ax, V3, V4)

    node(ax, G1, "", color=GRAY)
    node(ax, G2, "", color=GRAY)
    for xy, lab in [(V1, "1"), (V2, "2"), (V3, "3"), (V4, "4")]:
        node(ax, xy, lab)

    # Passenger groups: A by the base, B & C stacked beside node 2
    PGA = (3.5, 6.0)
    PGB = (10.7, 2.9)
    PGC = (10.7, 1.3)
    pgroup(ax, PGA, "PG A", "2 pax, direct-only", "1 \u2192 2")
    pgroup(ax, PGB, "PG B", "3 pax, stopover allowed", "2 \u2192 4")
    pgroup(ax, PGC, "PG C", "1 pax, stopover allowed", "2 \u2192 3")
    leader(ax, (PGA[0] + 0.2, PGA[1] - 0.5), V1)
    leader(ax, (PGB[0] - 1.35, PGB[1] - 0.2), (V2[0] + 0.2, V2[1] + 0.05))
    leader(ax, (PGC[0] - 1.35, PGC[1] + 0.2), (V2[0] + 0.2, V2[1] - 0.15))

    # Structural callouts
    callout(ax, (3.5, 3.45), "base vertiport", fs=11)
    callout(ax, (3.3, 1.95), "chosen end-vertiport", fs=11)
    leader(ax, (3.3, 2.15), (V4[0], V4[1] - 0.28))
    callout(ax, (3.3, 1.35),
            "region within 60 min drive of base  (valid end-vertiports)", fs=10.5)

    callout(ax, (-0.35, 3.5), "alternative", fs=10, ha="left")
    callout(ax, (-0.35, 3.2), "end-vertiport", fs=10, ha="left")
    leader(ax, (0.85, 3.35), G1)

    callout(ax, (7.55, 2.95), "shared leg", fs=10.5, ha="left", color=DARK)
    callout(ax, (6.55, 1.62), "(PG B + PG C)", fs=9.5, ha="left", color=DARK)
    leader(ax, (6.5, 2.5), (V2[0], 2.5))

    draw_legend(ax, 0.2, 0.5)

    ax.set_xlim(-1.4, 14.4)
    ax.set_ylim(-2.5, 7.5)
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