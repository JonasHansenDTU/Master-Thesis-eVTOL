from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go


def _safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def build_animation(snapshot_csv: Path, output_html: Path, title: str = "eVTOL Solution Animation") -> None:
    df = pd.read_csv(snapshot_csv)

    required = {
        "time",
        "evtol_id",
        "state",
        "node_from",
        "node_to",
        "op",
        "is_p",
        "is_o",
        "x",
        "y",
        "x_from",
        "y_from",
        "x_to",
        "y_to",
    }
    missing = sorted(required - set(df.columns))
    if missing:
        raise ValueError(f"Missing required columns in snapshot CSV: {missing}")

    # Normalize numeric columns
    for col in [
        "time",
        "evtol_id",
        "node_from",
        "node_to",
        "op",
        "is_p",
        "is_o",
        "x",
        "y",
        "x_from",
        "y_from",
        "x_to",
        "y_to",
    ]:
        df[col] = _safe_num(df[col])

    df["time"] = df["time"].astype("Int64")
    df["evtol_id"] = df["evtol_id"].astype("Int64")

    # Build static node layer from any known node endpoint coordinate
    node_points_a = df[["node_from", "x_from", "y_from"]].rename(
        columns={"node_from": "node", "x_from": "x", "y_from": "y"}
    )
    node_points_b = df[["node_to", "x_to", "y_to"]].rename(
        columns={"node_to": "node", "x_to": "x", "y_to": "y"}
    )
    nodes = (
        pd.concat([node_points_a, node_points_b], ignore_index=True)
        .dropna(subset=["node", "x", "y"])
        .drop_duplicates(subset=["node"])
        .sort_values("node")
    )

    times = sorted(df["time"].dropna().astype(int).unique().tolist())
    if not times:
        raise ValueError("No time values found in snapshot CSV.")

    evtols = sorted(df["evtol_id"].dropna().astype(int).unique().tolist())
    color_cycle = [
        "#1f77b4",
        "#d62728",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#8c564b",
        "#e377c2",
        "#17becf",
    ]
    evtol_color = {n: color_cycle[i % len(color_cycle)] for i, n in enumerate(evtols)}

    # Build animated point position for each (evtol, time):
    # - parked/inactive -> node coordinate
    # - flying -> linearly interpolated position from origin to destination
    #   across the contiguous flying time span for the same (evtol, op, arc).
    df = df.sort_values(["evtol_id", "time"]).reset_index(drop=True)
    df["anim_x"] = df["x"]
    df["anim_y"] = df["y"]

    # Parked points should sit exactly on the node.
    parked_mask = (df["is_p"].fillna(0.0) > 0.5) & df["x_from"].notna() & df["y_from"].notna()
    df.loc[parked_mask, "anim_x"] = df.loc[parked_mask, "x_from"]
    df.loc[parked_mask, "anim_y"] = df.loc[parked_mask, "y_from"]

    for n in evtols:
        dn = df[df["evtol_id"] == n].copy().sort_values("time")
        idx = dn.index.to_list()
        k = 0
        while k < len(idx):
            i0 = idx[k]
            row0 = df.loc[i0]
            flying0 = float(row0.get("is_o", 0.0) or 0.0) > 0.5
            valid_coords = pd.notna(row0["x_from"]) and pd.notna(row0["y_from"]) and pd.notna(row0["x_to"]) and pd.notna(row0["y_to"])
            if not (flying0 and valid_coords):
                k += 1
                continue

            # Build one contiguous flying segment with same arc/op and consecutive time.
            seg = [i0]
            cur_time = int(row0["time"])
            cur_op = row0["op"]
            cur_from = row0["node_from"]
            cur_to = row0["node_to"]

            j = k + 1
            while j < len(idx):
                ij = idx[j]
                rj = df.loc[ij]
                flyingj = float(rj.get("is_o", 0.0) or 0.0) > 0.5
                if not flyingj:
                    break
                same_trip = (
                    int(rj["time"]) == cur_time + 1
                    and rj["op"] == cur_op
                    and rj["node_from"] == cur_from
                    and rj["node_to"] == cur_to
                )
                if not same_trip:
                    break
                seg.append(ij)
                cur_time += 1
                j += 1

            x_from = float(row0["x_from"])
            y_from = float(row0["y_from"])
            x_to = float(row0["x_to"])
            y_to = float(row0["y_to"])

            if len(seg) == 1:
                df.loc[seg[0], "anim_x"] = x_from
                df.loc[seg[0], "anim_y"] = y_from
            else:
                for s_idx, row_idx in enumerate(seg):
                    frac = s_idx / (len(seg) - 1)
                    df.loc[row_idx, "anim_x"] = x_from + frac * (x_to - x_from)
                    df.loc[row_idx, "anim_y"] = y_from + frac * (y_to - y_from)

            k = j

    x_min = min(nodes["x"].min(), df["x_from"].min(), df["x_to"].min(), df["anim_x"].min())
    x_max = max(nodes["x"].max(), df["x_from"].max(), df["x_to"].max(), df["anim_x"].max())
    y_min = min(nodes["y"].min(), df["y_from"].min(), df["y_to"].min(), df["anim_y"].min())
    y_max = max(nodes["y"].max(), df["y_from"].max(), df["y_to"].max(), df["anim_y"].max())
    x_pad = (x_max - x_min) * 0.08 if x_max > x_min else 0.01
    y_pad = (y_max - y_min) * 0.08 if y_max > y_min else 0.01

    def traces_for_time(t: int) -> list[go.Scatter]:
        dft = df[df["time"] == t].copy()

        traces: list[go.Scatter] = []

        # Static nodes
        traces.append(
            go.Scatter(
                x=nodes["x"],
                y=nodes["y"],
                mode="markers+text",
                text=[str(int(v)) for v in nodes["node"]],
                textposition="top center",
                marker=dict(size=10, color="#444444"),
                name="Nodes",
                hovertemplate="Node %{text}<extra></extra>",
            )
        )

        # Per-eVTOL marker (always a dot), interpolated while flying.
        for n in evtols:
            row = dft[dft["evtol_id"] == n]
            if row.empty:
                continue
            r = row.iloc[0]
            color = evtol_color[n]

            is_flying = float(r.get("is_o", 0.0) or 0.0) > 0.5
            is_parked = float(r.get("is_p", 0.0) or 0.0) > 0.5

            # Vehicle marker at current representative position
            if pd.notna(r["anim_x"]) and pd.notna(r["anim_y"]):
                state = "flying" if is_flying else ("parked" if is_parked else "inactive")
                traces.append(
                    go.Scatter(
                        x=[r["anim_x"]],
                        y=[r["anim_y"]],
                        mode="markers+text",
                        text=[f"{n}"],
                        textposition="top center",
                        marker=dict(size=14, color=color, symbol="circle", line=dict(width=1, color="#222")),
                        name=f"eVTOL {n}",
                        hovertemplate=(
                            f"eVTOL {n}<br>state={state}"
                            + f"<br>time={t}"
                            + f"<br>node_from={int(r['node_from']) if pd.notna(r['node_from']) else '-'}"
                            + f"<br>node_to={int(r['node_to']) if pd.notna(r['node_to']) else '-'}"
                            + f"<br>op={int(r['op']) if pd.notna(r['op']) else '-'}"
                            + f"<br>is_p={float(r['is_p']):.2f}"
                            + f"<br>is_o={float(r['is_o']):.2f}"
                            + "<extra></extra>"
                        ),
                        showlegend=False,
                    )
                )

        return traces

    initial_traces = traces_for_time(times[0])
    frames = [go.Frame(data=traces_for_time(t), name=str(t)) for t in times]

    fig = go.Figure(data=initial_traces, frames=frames)
    fig.update_layout(
        title=f"{title} | t = {times[0]}",
        template="plotly_white",
        xaxis=dict(title="Longitude", range=[x_min - x_pad, x_max + x_pad], scaleanchor="y", scaleratio=1),
        yaxis=dict(title="Latitude", range=[y_min - y_pad, y_max + y_pad]),
        width=1100,
        height=750,
        margin=dict(l=40, r=20, t=70, b=40),
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                x=0.02,
                y=1.12,
                xanchor="left",
                yanchor="top",
                buttons=[
                    dict(label="Play", method="animate", args=[None, {"frame": {"duration": 250, "redraw": True}, "fromcurrent": True}]),
                    dict(label="Pause", method="animate", args=[[None], {"frame": {"duration": 0, "redraw": False}, "mode": "immediate"}]),
                ],
            )
        ],
        sliders=[
            dict(
                active=0,
                x=0.12,
                y=1.12,
                len=0.84,
                xanchor="left",
                yanchor="top",
                currentvalue={"prefix": "time = "},
                steps=[
                    {
                        "label": str(t),
                        "method": "animate",
                        "args": [[str(t)], {"frame": {"duration": 0, "redraw": True}, "mode": "immediate"}],
                    }
                    for t in times
                ],
            )
        ],
    )

    # Dynamic title per frame via layout update in slider steps is optional; keep static + slider current value.
    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")
    print(f"Animation written: {output_html}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Plotly slider animation from eVTOL solution snapshots.")
    parser.add_argument(
        "--input",
        type=Path,
        default=Path(__file__).with_name("solution_snapshots.csv"),
        help="Path to snapshot CSV exported from Model.jl",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path(__file__).with_name("solution_animation.html"),
        help="Output HTML animation path",
    )
    parser.add_argument("--title", type=str, default="eVTOL Solution Animation", help="Figure title")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    build_animation(args.input, args.output, title=args.title)
