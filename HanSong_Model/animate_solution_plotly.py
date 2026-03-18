from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots

try:
    import folium
    from folium.plugins import TimestampedGeoJson
except ImportError:  # pragma: no cover
    folium = None
    TimestampedGeoJson = None


def _safe_num(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _nearest_valid(values: list[float], idx: int, direction: int) -> float | None:
    i = idx
    while 0 <= i < len(values):
        v = values[i]
        if not math.isnan(v):
            return v
        i += direction
    return None


def _prepare_snapshot_data(snapshot_csv: Path, infrastructure_csv: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame, list[int], list[int]]:
    """Read and normalize snapshot CSV into a dataframe ready for visualization.

    Returns:
        df: Full snapshot dataframe with interpolated coordinates (`anim_x`, `anim_y`) and battery smoothing.
        nodes: Unique nodes extracted from `node_from/node_to` endpoints.
        times: Sorted list of unique integer time steps.
        evtols: Sorted list of unique eVTOL identifiers.
    """

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

    # Optional enhancement columns from newer snapshot exports.
    if "battery_level" in df.columns:
        df["battery_level"] = _safe_num(df["battery_level"])
    else:
        df["battery_level"] = float("nan")

    if "onboard_passenger_count" in df.columns:
        df["onboard_passenger_count"] = _safe_num(df["onboard_passenger_count"])
    else:
        df["onboard_passenger_count"] = float("nan")

    if "onboard_groups" not in df.columns:
        df["onboard_groups"] = ""

    # Build a smoothed battery series at integer times under a linear assumption:
    # - flying segments: linear decrease
    # - ground segments (parked/inactive): linear increase or flat
    # This removes staircase artifacts from operation-level battery reporting.
    df["battery_smooth"] = df["battery_level"]

    df["time"] = df["time"].astype("Int64")
    df["evtol_id"] = df["evtol_id"].astype("Int64")

    # Build static node layer
    if infrastructure_csv is not None and infrastructure_csv.exists():
        infra_df = pd.read_csv(infrastructure_csv)
        # Assume columns: node, lat, lon
        infra_df = infra_df.rename(columns={"lon": "x", "lat": "y"})
        nodes = infra_df[["node", "x", "y"]].dropna().drop_duplicates(subset=["node"]).sort_values("node")
    else:
        # Build from snapshot data
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

        # 1) Existing position interpolation prep.
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

        # 2) Battery smoothing across mode segments.
        batt_raw = [float(df.loc[i, "battery_level"]) if pd.notna(df.loc[i, "battery_level"]) else float("nan") for i in idx]

        # Fill hard NaNs by nearest values so interpolation has anchors.
        for p in range(len(batt_raw)):
            if math.isnan(batt_raw[p]):
                left = _nearest_valid(batt_raw, p - 1, -1)
                right = _nearest_valid(batt_raw, p + 1, +1)
                if left is not None and right is not None:
                    batt_raw[p] = 0.5 * (left + right)
                elif left is not None:
                    batt_raw[p] = left
                elif right is not None:
                    batt_raw[p] = right

        mode = []
        for i in idx:
            is_o = float(df.loc[i, "is_o"] if pd.notna(df.loc[i, "is_o"]) else 0.0)
            mode.append("flying" if is_o > 0.5 else "ground")

        s = 0
        while s < len(idx):
            e = s
            while e + 1 < len(idx) and mode[e + 1] == mode[s]:
                e += 1

            # Segment anchors from neighboring points when possible.
            prev_b = batt_raw[s - 1] if s - 1 >= 0 else batt_raw[s]
            next_b = batt_raw[e + 1] if e + 1 < len(idx) else batt_raw[e]

            if mode[s] == "flying":
                b0 = prev_b
                b1 = next_b
                # Flying should not increase under linear usage assumption.
                if b1 > b0:
                    b1 = b0
            else:
                b0 = prev_b
                b1 = next_b
                # Ground charging should not decrease under linear assumption.
                if b1 < b0:
                    b1 = b0

            seg_len = e - s + 1
            for k_seg in range(seg_len):
                frac = 0.0 if seg_len == 1 else k_seg / (seg_len - 1)
                b = (1.0 - frac) * b0 + frac * b1
                df.loc[idx[s + k_seg], "battery_smooth"] = b

            s = e + 1

    return df, nodes, times, evtols


def build_animation(
    snapshot_csv: Path,
    output_html: Path,
    title: str = "eVTOL Solution Animation",
    fps: float = 12.0,
    subframes: int = 1,
    infrastructure_csv: Path | None = None,
) -> None:
    if fps <= 0:
        raise ValueError("fps must be > 0")
    if subframes < 1:
        raise ValueError("subframes must be >= 1")

    # Keep the same total runtime per original time step.
    # Example: subframes=4 -> 4x more frames, each with 1/4 duration.
    base_frame_duration_ms = max(1, int(round(1000.0 / fps)))
    frame_duration_ms = max(1, int(round(base_frame_duration_ms / subframes)))

    df, nodes, times, evtols = _prepare_snapshot_data(snapshot_csv, infrastructure_csv)

    color_cycle = [
        "#1f77b4",
        "#cb2020",
        "#2ca02c",
        "#ff7f0e",
        "#9467bd",
        "#ffe100",
        "#e377c2",
        "#17becf",
    ]
    evtol_color = {n: color_cycle[i % len(color_cycle)] for i, n in enumerate(evtols)}

    x_min = min(nodes["x"].min(), df["x_from"].min(), df["x_to"].min(), df["anim_x"].min())
    x_max = max(nodes["x"].max(), df["x_from"].max(), df["x_to"].max(), df["anim_x"].max())
    y_min = min(nodes["y"].min(), df["y_from"].min(), df["y_to"].min(), df["anim_y"].min())
    y_max = max(nodes["y"].max(), df["y_from"].max(), df["y_to"].max(), df["anim_y"].max())
    x_pad = (x_max - x_min) * 0.08 if x_max > x_min else 0.01
    y_pad = (y_max - y_min) * 0.08 if y_max > y_min else 0.01
    battery_text_offset = (y_max - y_min) * 0.015 if y_max > y_min else 0.005

    # Fast lookup for interpolation.
    row_lookup: dict[tuple[int, int], pd.Series] = {
        (int(r["evtol_id"]), int(r["time"])): r for _, r in df.iterrows()
    }

    def _interp_scalar(a: float, b: float, alpha: float) -> float:
        return (1.0 - alpha) * a + alpha * b

    def _snapshot_for_evtol(n: int, tau: float):
        t_floor = math.floor(tau)
        t_ceil = math.ceil(tau)
        alpha = 0.0 if t_ceil == t_floor else (tau - t_floor) / (t_ceil - t_floor)
        r0 = row_lookup.get((n, t_floor))
        r1 = row_lookup.get((n, t_ceil))
        if r0 is None and r1 is None:
            return None
        if r0 is None:
            r0 = r1
        if r1 is None:
            r1 = r0

        x0 = float(r0["anim_x"])
        y0 = float(r0["anim_y"])
        x1 = float(r1["anim_x"])
        y1 = float(r1["anim_y"])
        x = _interp_scalar(x0, x1, alpha)
        y = _interp_scalar(y0, y1, alpha)

        is_p = _interp_scalar(float(r0.get("is_p", 0.0) or 0.0), float(r1.get("is_p", 0.0) or 0.0), alpha)
        is_o = _interp_scalar(float(r0.get("is_o", 0.0) or 0.0), float(r1.get("is_o", 0.0) or 0.0), alpha)
        b0 = float(r0["battery_smooth"]) if pd.notna(r0.get("battery_smooth")) else float("nan")
        b1 = float(r1["battery_smooth"]) if pd.notna(r1.get("battery_smooth")) else float("nan")
        battery = _interp_scalar(b0, b1, alpha) if (not math.isnan(b0) and not math.isnan(b1)) else (b0 if not math.isnan(b0) else b1)

        ref_row = r0 if alpha < 0.5 else r1
        pax = int(ref_row["onboard_passenger_count"]) if pd.notna(ref_row.get("onboard_passenger_count")) else 0
        groups = str(ref_row.get("onboard_groups", "") or "")
        state = "flying" if is_o > 0.5 else ("parked" if is_p > 0.5 else "inactive")

        return {
            "x": x,
            "y": y,
            "is_p": is_p,
            "is_o": is_o,
            "battery": battery,
            "state": state,
            "pax": pax,
            "groups": groups,
            "ref_row": ref_row,
        }

    def map_traces_for_time(tau: float) -> list[go.Scatter]:
        traces: list[go.Scatter] = []

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

        for n in evtols:
            snap = _snapshot_for_evtol(n, tau)
            if snap is None:
                traces.append(
                    go.Scatter(
                        x=[None],
                        y=[None],
                        mode="markers",
                        marker=dict(size=14, color=evtol_color[n], symbol="circle", line=dict(width=1, color="#222")),
                        showlegend=False,
                        name=f"eVTOL {n}",
                    )
                )
                traces.append(go.Scatter(x=[None], y=[None], mode="text", showlegend=False, hoverinfo="skip"))
                continue

            color = evtol_color[n]
            ref_row = snap["ref_row"]
            battery = snap["battery"]

            traces.append(
                go.Scatter(
                    x=[snap["x"]],
                    y=[snap["y"]],
                    mode="markers",
                    marker=dict(size=14, color=color, symbol="circle", line=dict(width=1, color="#222")),
                    name=f"eVTOL {n}",
                    hovertemplate=(
                        f"eVTOL {n}<br>state={snap['state']}"
                        + f"<br>time={tau:.2f}"
                        + f"<br>node_from={int(ref_row['node_from']) if pd.notna(ref_row['node_from']) else '-'}"
                        + f"<br>node_to={int(ref_row['node_to']) if pd.notna(ref_row['node_to']) else '-'}"
                        + f"<br>op={int(ref_row['op']) if pd.notna(ref_row['op']) else '-'}"
                        + f"<br>is_p={snap['is_p']:.2f}"
                        + f"<br>is_o={snap['is_o']:.2f}"
                        + f"<br>pax={snap['pax']}"
                        + (f"<br>groups={snap['groups']}" if snap['groups'] else "")
                        + (f"<br>battery={battery:.1f}%" if not math.isnan(battery) else "")
                        + "<extra></extra>"
                    ),
                    showlegend=False,
                )
            )

            if not math.isnan(battery):
                traces.append(
                    go.Scatter(
                        x=[snap["x"]],
                        y=[snap["y"] + battery_text_offset],
                        mode="text",
                        text=[f"{battery:.0f}%"],
                        textposition="middle center",
                        textfont=dict(size=10, color=color),
                        hoverinfo="skip",
                        showlegend=False,
                    )
                )
            else:
                traces.append(go.Scatter(x=[None], y=[None], mode="text", showlegend=False, hoverinfo="skip"))

        return traces

    def table_trace_for_time(tau: float) -> go.Table:
        rows = []
        for n in evtols:
            snap = _snapshot_for_evtol(n, tau)
            if snap is None:
                rows.append((str(n), "-", "-", "-", "-"))
                continue
            batt_txt = "-" if math.isnan(snap["battery"]) else f"{snap['battery']:.1f}%"
            groups_txt = snap["groups"] if snap["groups"] and snap["groups"].lower() != "nan" else "-"
            rows.append((str(n), snap["state"], str(snap["pax"]), batt_txt, groups_txt))

        col_evtol = [r[0] for r in rows]
        col_state = [r[1] for r in rows]
        col_pax = [r[2] for r in rows]
        col_batt = [r[3] for r in rows]
        col_groups = [r[4] for r in rows]

        return go.Table(
            header=dict(
                values=["eVTOL", "State", "Pax", "Battery", "Groups"],
                fill_color="#e9eef5",
                align="left",
                font=dict(size=12),
            ),
            cells=dict(
                values=[col_evtol, col_state, col_pax, col_batt, col_groups],
                align="left",
                font=dict(size=11),
                fill_color="#ffffff",
                height=26,
            ),
            columnwidth=[60, 75, 50, 80, 260],
        )

    frame_times: list[float] = []
    if len(times) == 1:
        frame_times = [float(times[0])]
    else:
        for i in range(len(times) - 1):
            t0 = float(times[i])
            t1 = float(times[i + 1])
            for s in range(subframes):
                frac = s / subframes
                frame_times.append(t0 + frac * (t1 - t0))
        frame_times.append(float(times[-1]))

    def _fmt_time_label(v: float) -> str:
        if abs(v - round(v)) < 1e-9:
            return str(int(round(v)))
        return f"{v:.2f}".rstrip("0").rstrip(".")

    initial_map_traces = map_traces_for_time(frame_times[0])
    initial_table_trace = table_trace_for_time(frame_times[0])
    initial_data = initial_map_traces + [initial_table_trace]

    total_traces = len(initial_data)
    frames = [
        go.Frame(
            data=map_traces_for_time(tau) + [table_trace_for_time(tau)],
            traces=list(range(total_traces)),
            name=_fmt_time_label(tau),
        )
        for tau in frame_times
    ]

    play_args = {
        "frame": {"duration": frame_duration_ms, "redraw": True},
        "transition": {"duration": 0},
        "mode": "immediate",
        "fromcurrent": True,
    }
    jump_args = {
        "frame": {"duration": 0, "redraw": True},
        "transition": {"duration": 0},
        "mode": "immediate",
    }

    fig = make_subplots(
        rows=1,
        cols=2,
        specs=[[{"type": "xy"}, {"type": "table"}]],
        column_widths=[0.72, 0.28],
        horizontal_spacing=0.03,
    )
    for tr in initial_data[:-1]:
        fig.add_trace(tr, row=1, col=1)
    fig.add_trace(initial_data[-1], row=1, col=2)
    fig.frames = frames
    fig.update_layout(
        title=f"{title} | t = {_fmt_time_label(frame_times[0])}",
        title_x=0.28,
        template="plotly_white",
        xaxis=dict(title="Longitude", range=[x_min - x_pad, x_max + x_pad], scaleanchor="y", scaleratio=1),
        yaxis=dict(title="Latitude", range=[y_min - y_pad, y_max + y_pad]),
        showlegend=False,
        width=1420,
        height=860,
        margin=dict(l=40, r=30, t=60, b=140),
        annotations=[
            dict(
                text="eVTOL Status",
                x=0.87,
                y=1.03,
                xref="paper",
                yref="paper",
                showarrow=False,
                font=dict(size=14),
            )
        ],
        updatemenus=[
            dict(
                type="buttons",
                showactive=False,
                x=0.03,
                y=-0.095,
                xanchor="left",
                yanchor="bottom",
                buttons=[
                    dict(
                        label="Play",
                        method="animate",
                        args=[None, play_args],
                    ),
                    dict(label="Pause", method="animate", args=[[None], jump_args]),
                ],
            )
        ],
        sliders=[
            dict(
                active=0,
                x=0.16,
                y=-0.1,
                len=0.52,
                xanchor="left",
                yanchor="bottom",
                pad={"t": 10, "b": 0},
                currentvalue={"prefix": "time = "},
                steps=[
                    {
                        "label": _fmt_time_label(tau),
                        "method": "animate",
                        "args": [[_fmt_time_label(tau)], jump_args],
                    }
                    for tau in frame_times
                ],
            )
        ],
    )

    fig.update_xaxes(title_text="Longitude", range=[x_min - x_pad, x_max + x_pad], scaleanchor="y", scaleratio=1, row=1, col=1)
    fig.update_yaxes(title_text="Latitude", range=[y_min - y_pad, y_max + y_pad], row=1, col=1)

    output_html.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_html), include_plotlyjs="cdn")
    print(f"Animation written: {output_html}")


def build_folium_map(
    snapshot_csv: Path,
    output_html: Path,
    title: str = "eVTOL Solution Map",
    start_time: str = "2020-01-01T00:00:00Z",
    time_unit: str = "s",
    infrastructure_csv: Path | None = None,
) -> None:
    """Generate a Folium map with a time slider showing eVTOL positions."""

    if folium is None or TimestampedGeoJson is None:
        raise ImportError(
            "folium is required for the folium backend. Install it with `pip install folium`."
        )

    df, nodes, _, _ = _prepare_snapshot_data(snapshot_csv, infrastructure_csv)

    # Prefer already-processed interpolated coords when available.
    x_col = "anim_x" if "anim_x" in df.columns else "x"
    y_col = "anim_y" if "anim_y" in df.columns else "y"

    # Assign a consistent color per eVTOL.
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

    # Time -> ISO string for TimestampedGeoJson.
    origin = pd.Timestamp(start_time)
    if origin.tzinfo is not None:
        origin = origin.tz_convert(None)
    df["timestamp"] = pd.to_datetime(df["time"], unit=time_unit, origin=origin)
    df["timestamp_str"] = df["timestamp"].dt.strftime("%Y-%m-%dT%H:%M:%SZ")

    # Build GeoJSON features for each snapshot.
    features = []
    status_by_time: dict[str, dict[str, dict[str, str]]] = {}

    for _, r in df.iterrows():
        lon = float(r[x_col])
        lat = float(r[y_col])
        if math.isnan(lon) or math.isnan(lat):
            continue

        time_str = r["timestamp_str"]
        evtol_id = str(int(r["evtol_id"])) if pd.notna(r["evtol_id"]) else None
        if evtol_id is not None:
            status_by_time.setdefault(time_str, {})[evtol_id] = {
                "state": str(r["state"]),
                "pax": str(int(r["onboard_passenger_count"])) if pd.notna(r.get("onboard_passenger_count")) else "-",
                "battery": f"{float(r['battery_smooth']):.1f}%" if pd.notna(r.get("battery_smooth")) else "-",
                "groups": str(r.get("onboard_groups") or "-"),
            }

        popup = (
            f"eVTOL {evtol_id if evtol_id is not None else '?'}<br>"
            f"time={time_str}<br>"
            f"state={r['state']}<br>"
            f"node_from={int(r['node_from']) if pd.notna(r['node_from']) else '-'}<br>"
            f"node_to={int(r['node_to']) if pd.notna(r['node_to']) else '-'}"
        )

        fill_color = evtol_color.get(int(r["evtol_id"]) if pd.notna(r["evtol_id"]) else None, "#1f77b4")
        properties = {
            "time": time_str,
            "popup": popup,
            "icon": "circle",
            "iconstyle": {
                "fillColor": fill_color,
                "fillOpacity": 0.8,
                "stroke": "true",
                "radius": 6,
            },
        }

        features.append(
            {
                "type": "Feature",
                "geometry": {"type": "Point", "coordinates": [lon, lat]},
                "properties": properties,
            }
        )

    geojson = {"type": "FeatureCollection", "features": features}

    if not df.empty:
        center_lat = float(df[y_col].dropna().mean())
        center_lon = float(df[x_col].dropna().mean())
    else:
        center_lat = 0.0
        center_lon = 0.0

    m = folium.Map(location=[center_lat, center_lon], zoom_start=11, tiles="OpenStreetMap")
    folium.TileLayer("cartodbpositron").add_to(m)

    for _, r in nodes.iterrows():
        folium.CircleMarker(
            location=[float(r["y"]), float(r["x"])],
            radius=5,
            color="#222222",
            fill=True,
            fill_color="#888888",
            fill_opacity=0.7,
            popup=f"Node {int(r['node'])}",
        ).add_to(m)

    # Add an info table overlay.
    evtol_ids = sorted({int(k) for t in status_by_time.values() for k in t.keys()})
    status_table_html = """
    <div id="evtol-status" style="position:absolute; top:10px; left:10px; z-index:9999; background:rgba(255,255,255,0.85); padding:10px; border-radius:8px; max-height:240px; overflow:auto; font-size:12px;">
      <div style="font-weight:bold; margin-bottom:6px;">eVTOL status</div>
      <table style="border-collapse:collapse; width:100%;">
        <thead>
          <tr>
            <th style="text-align:left; padding:2px 4px;">eVTOL</th>
            <th style="text-align:left; padding:2px 4px;">State</th>
            <th style="text-align:right; padding:2px 4px;">Pax</th>
            <th style="text-align:right; padding:2px 4px;">Battery</th>
            <th style="text-align:left; padding:2px 4px;">Groups</th>
          </tr>
        </thead>
        <tbody id="evtol-status-table-body">
        </tbody>
      </table>
    </div>
    """
    m.get_root().html.add_child(folium.Element(status_table_html))

    TimestampedGeoJson(
        geojson,
        period="PT1S",
        add_last_point=True,
        duration="PT2S",
        auto_play=False,
        loop=False,
        max_speed=1,
        loop_button=True,
        date_options="YYYY/MM/DD HH:mm:ss",
        time_slider_drag_update=True,
    ).add_to(m)

    # Inject JS to update table as the timeline changes.
    import json

    js = """
    <script>
    const evtolStatus = %s;
    const mapName = "%s";

    function updateEVTOLTable(timestamp) {
      const body = document.getElementById('evtol-status-table-body');
      if (!body) return;
      const status = evtolStatus[timestamp] || {};
      const rows = [];
      const ids = %s;
      for (const id of ids) {
        const s = status[id] || {state: '-', pax: '-', battery: '-', groups: '-'};
        rows.push(
          `<tr>` +
            `<td style="padding:2px 4px;">${id}</td>` +
            `<td style="padding:2px 4px;">${s.state}</td>` +
            `<td style="padding:2px 4px; text-align:right;">${s.pax}</td>` +
            `<td style="padding:2px 4px; text-align:right;">${s.battery}</td>` +
            `<td style="padding:2px 4px;">${s.groups}</td>` +
          `</tr>`
        );
      }
      body.innerHTML = rows.join('');
    }

    function toTimestampKey(v) {
      // Normalize to the exact keys we generate in Python: YYYY-MM-DDTHH:MM:SSZ
      const date = typeof v === 'string' ? new Date(v) : new Date(v);
      if (Number.isNaN(date.getTime())) {
        return v;
      }
      return date.toISOString().slice(0,19) + 'Z';
    }

    window.addEventListener('load', () => {
      const mapObj = window[mapName];
      if (!mapObj || !mapObj.timeDimension) {
        console.warn('Timedimension not found for map:', mapName);
        return;
      }

      mapObj.timeDimension.on('timeload', (e) => {
        updateEVTOLTable(toTimestampKey(e.time));
      });

      // Initialize table immediately.
      const initialTime = Object.keys(evtolStatus)[0];
      if (initialTime) updateEVTOLTable(initialTime);
    });
    </script>
    """ % (
        json.dumps(status_by_time),
        m.get_name(),
        json.dumps(evtol_ids),
    )

    m.get_root().html.add_child(folium.Element(js))

    output_html.parent.mkdir(parents=True, exist_ok=True)
    m.save(str(output_html))
    print(f"Folium map written: {output_html}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build an animated view of an eVTOL solution (Plotly or Folium).")
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
        help="Output HTML animation/map path",
    )
    parser.add_argument(
        "--backend",
        type=str,
        choices=["plotly", "folium"],
        default="folium",
        help="Visualization backend to use (plotly or folium).",
    )
    parser.add_argument("--title", type=str, default="eVTOL Solution Visualization", help="Figure title")
    parser.add_argument(
        "--fps",
        type=float,
        default=14.0,
        help="Base playback FPS per original integer time-step (default: 14)",
    )
    parser.add_argument(
        "--subframes",
        type=int,
        default=5,
        help="Interpolated frames inserted between consecutive time steps while keeping total runtime per step unchanged",
    )
    parser.add_argument(
        "--start-time",
        type=str,
        default="2020-01-01T00:00:00Z",
        help="Base timestamp used for folium timeline (time=0 corresponds to this).",
    )
    parser.add_argument(
        "--infrastructure",
        type=Path,
        default=Path(__file__).with_name("infrastructure.csv"),
        help="Path to infrastructure CSV with columns: node, lat, lon. Shows all vertiports even if not used in solution.",
    )
    parser.add_argument(
        "--time-unit",
        type=str,
        default="s",
        choices=["s", "m", "h", "d"],
        help="Unit of the numeric 'time' column when generating folium timestamps (s=seconds, m=minutes, h=hours, d=days).",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.backend == "plotly":
        build_animation(args.input, args.output, title=args.title, fps=args.fps, subframes=args.subframes, infrastructure_csv=args.infrastructure)
    else:
        build_folium_map(
            args.input,
            args.output,
            title=args.title,
            start_time=args.start_time,
            time_unit=args.time_unit,
            infrastructure_csv=args.infrastructure,
        )
