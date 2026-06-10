import pandas as pd
import streamlit as st
import plotly.express as px

def run_dashboard(df: pd.DataFrame):
    # --- Page Config & High-Contrast CSS Theme ---
    st.set_page_config(page_title="eVTOL Optimization Control Room", layout="wide", page_icon="🛸")
    
    # Custom CSS for structure, deep text-to-background contrast, and modern layout cards
    st.markdown("""
        <style>
            /* Reset canvas padding */
            .block-container { padding-top: 1.5rem; padding-bottom: 2rem; background-color: #FAFAFA; }
            
            /* High contrast typography (Dark Slate on Off-White) */
            h1, h2, h3 { color: #0F172A !important; font-weight: 700 !important; }
            p, span, label { color: #334155 !important; font-weight: 500; }
            
            /* Custom KPI Card Styling with high-contrast borders */
            .kpi-card {
                background-color: #FFFFFF;
                padding: 20px;
                border-radius: 12px;
                border: 2px solid #E2E8F0;
                box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.05);
                margin-bottom: 1rem;
            }
            .kpi-title { font-size: 0.85rem; text-transform: uppercase; letter-spacing: 0.05em; color: #64748B; margin-bottom: 6px; }
            .kpi-value { font-size: 1.8rem; font-weight: 800; color: #0F172A; }
            
            /* Section container */
            .dashboard-section {
                background-color: #FFFFFF;
                padding: 24px;
                border-radius: 16px;
                border: 1px solid #E2E8F0;
                margin-bottom: 1.5rem;
            }
        </style>
    """, unsafe_allow_html=True)

    # --- Header ---
    st.title("🛸 eVTOL Fleet Optimization Control Room")
    st.markdown("Assess solver performance, vehicle bottlenecks, and passenger service levels across execution history.")
    st.markdown("---")

    # --- Sidebar Filters ---
    st.sidebar.header("🎛️ Global Controls")
    
    run_min, run_max = int(df["run"].min()), int(df["run"].max())
    run_range = st.sidebar.slider(
        "Filter Run Execution Window",
        run_min, run_max, (run_min, run_max)
    )
    
    # Filter and sort data
    filtered = df[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])].sort_values("run")

    # --- SECTION 1: EXECUTIVE SCORECARD ---
    st.subheader("📊 Executive Summary Metrics (Averages)")
    
    # Compute metrics cleanly upfront
    max_obj = filtered['best_obj'].max() if 'best_obj' in filtered.columns else 0.0
    avg_util = (filtered['aircraft_utilization'].mean() * 100) if 'aircraft_utilization' in filtered.columns else 0.0
    avg_demand = filtered['passenger_demand_served'].mean() if 'passenger_demand_served' in filtered.columns else 0.0
    avg_wait = filtered['average_passenger_waiting_time'].mean() if 'average_passenger_waiting_time' in filtered.columns else 0.0
    avg_deadhead = filtered['deadhead_time'].mean() if 'deadhead_time' in filtered.columns else 0.0

    # Grid of HTML-rendered High Contrast KPI Cards
    kpi_cols = st.columns(5)
    
    metrics_data = [
        ("🏆 Max Best Objective", f"{max_obj:,.2f}", 0),
        ("⚡ Avg Utilization", f"{avg_util:.1f}%", 1),
        ("👥 Avg Demand Served", f"{avg_demand:,.1f}", 2),
        ("⏳ Avg Waiting Time", f"{avg_wait:.1f} min", 3),
        ("✈️ Avg Deadhead Time", f"{avg_deadhead:.1f} min", 4)
    ]
    
    for title, val, col_idx in metrics_data:
        kpi_cols[col_idx].markdown(f"""
            <div class="kpi-card">
                <div class="kpi-title">{title}</div>
                <div class="kpi-value">{val}</div>
            </div>
        """, unsafe_allow_html=True)

    # --- SECTION 2: MACRO ANALYSIS (SIDE-BY-SIDE) ---
    st.markdown('<div class="dashboard-section">', unsafe_allow_html=True)
    mid_left, mid_right = st.columns([3, 2])
    
    target_kpis = ["best_obj", "aircraft_utilization", "deadhead_time", "passenger_demand_served", "average_passenger_waiting_time"]
    available_kpis = [k for k in target_kpis if k in df.columns]

    with mid_left:
        st.markdown("### 📊 Distribution of Best Objective Values")
        if "best_obj" in filtered.columns:
            # Replaced px.line with px.box for distributional overview
            fig_main = px.box(
                filtered, 
                y="best_obj",
                points="all", # Shows individual run dots next to the box for transparency
                template="plotly_white", 
                color_discrete_sequence=["#2563EB"] # Bold blue theme
            )
            fig_main.update_layout(
                height=350, 
                margin=dict(l=10, r=10, t=20, b=10),
                # Enforce clean high-contrast light theme attributes
                paper_bgcolor="#FFFFFF",
                plot_bgcolor="#FFFFFF",
                font=dict(color="#0F172A"),
                xaxis=dict(
                    showticklabels=False, # Hides the default empty X-axis label for a single box
                    gridcolor="#E2E8F0"
                ),
                yaxis=dict(
                    gridcolor="#E2E8F0", 
                    title="Objective Value Value",
                    title_font=dict(color="#0F172A"), 
                    tickfont=dict(color="#334155")
                )
            )
            st.plotly_chart(fig_main, use_container_width=True)

    with mid_right:
        st.markdown("### 🎯 KPI Trade-Off & Correlation")
        if len(available_kpis) > 1:
            corr = filtered[available_kpis].corr()
            
            # Using an explicitly distinct light-to-dark blue gradient map for high contrast
            fig_corr = px.imshow(
                corr, text_auto=".2f", aspect="auto",
                template="plotly_white",
                color_continuous_scale=[[0, "#EFF6FF"], [0.5, "#3B82F6"], [1, "#1E3A8A"]]
            )
            fig_corr.update_layout(
                height=350, 
                margin=dict(l=10, r=10, t=20, b=10),
                # Force high-contrast light theme attributes
                paper_bgcolor="#FFFFFF",
                plot_bgcolor="#FFFFFF",
                font=dict(color="#0F172A"),
                coloraxis_showscale=False,
                xaxis=dict(
                    title_font=dict(color="#0F172A"), 
                    tickfont=dict(color="#334155")
                ),
                yaxis=dict(
                    title_font=dict(color="#0F172A"), 
                    tickfont=dict(color="#334155")
                )
            )
            st.plotly_chart(fig_corr, use_container_width=True)

    # --- SECTION 3: EFFICIENCY & GRIDS DEEP DIVE ---
    st.subheader("⚙️ Operational Efficiency Deep-Dive")
    
    bot_col1, bot_col2 = st.columns(2)
    
    with bot_col1:
        # Using native containers with borders for high contrast and flawless rendering
        with st.container(border=True):
            st.markdown("### ✈️ Fleet Utilization vs. Empty Transit")
            
            # Dual visualization tracking operations
            if "aircraft_utilization" in filtered.columns and "deadhead_time" in filtered.columns:
                fig_fleet = px.scatter(
                    filtered, 
                    x="aircraft_utilization", 
                    y="deadhead_time", 
                    color="run", 
                    size="passenger_demand_served" if "passenger_demand_served" in filtered.columns else None,
                    labels={"aircraft_utilization": "Aircraft Utilization (Ratio)", "deadhead_time": "Deadhead Time (Min)", "run": "Run ID"},
                    template="plotly_white", 
                    color_continuous_scale="Viridis"
                )
                fig_fleet.update_layout(
                    height=350, 
                    margin=dict(l=10, r=10, t=20, b=10),
                    # Force high-contrast text and clean light canvas background
                    paper_bgcolor="#FFFFFF",
                    plot_bgcolor="#FFFFFF",
                    font=dict(color="#0F172A"), # Deep slate text
                    xaxis=dict(
                        gridcolor="#E2E8F0", 
                        title_font=dict(color="#0F172A"), 
                        tickfont=dict(color="#334155")
                    ),
                    yaxis=dict(
                        gridcolor="#E2E8F0", 
                        title_font=dict(color="#0F172A"), 
                        tickfont=dict(color="#334155")
                    )
                )
                st.plotly_chart(fig_fleet, use_container_width=True)
            else:
                st.info("Missing 'aircraft_utilization' or 'deadhead_time' columns to render plot.")

    with bot_col2:
        with st.container(border=True):
            st.markdown("### 👥 Passenger Experience Metrics")
            
            # Line/Area chart for passenger wait times over execution runs
            if "average_passenger_waiting_time" in filtered.columns:
                fig_pass = px.area(
                    filtered, 
                    x="run", 
                    y="average_passenger_waiting_time",
                    labels={"run": "Run Number", "average_passenger_waiting_time": "Avg Wait Time (Min)"},
                    template="plotly_white", 
                    color_discrete_sequence=["#F59E0B"] # Sharp amber color
                )
                fig_pass.update_layout(
                    hovermode="x unified", 
                    height=350, 
                    margin=dict(l=10, r=10, t=20, b=10),
                    # Force high-contrast text and clean light canvas background
                    paper_bgcolor="#FFFFFF",
                    plot_bgcolor="#FFFFFF",
                    font=dict(color="#0F172A"),
                    xaxis=dict(
                        gridcolor="#E2E8F0", 
                        title="Run Number",
                        title_font=dict(color="#0F172A"), 
                        tickfont=dict(color="#334155")
                    ),
                    yaxis=dict(
                        gridcolor="#E2E8F0", 
                        title="Avg Wait Time (Min)",
                        title_font=dict(color="#0F172A"), 
                        tickfont=dict(color="#334155")
                    )
                )
                st.plotly_chart(fig_pass, use_container_width=True)
            else:
                st.info("Missing 'average_passenger_waiting_time' column to render plot.")