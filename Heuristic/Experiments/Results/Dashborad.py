import pandas as pd
import streamlit as st
import plotly.express as px

def run_dashboard(df: pd.DataFrame):
    st.set_page_config(page_title="Run Dashboard", layout="wide")

    st.title("Optimization Run Dashboard")

    # --- Sidebar controls ---
    metric_options = [c for c in df.columns if c != "run"]
    selected_metrics = st.sidebar.multiselect(
        "Select metrics to visualize",
        metric_options,
        default=metric_options[:4]
    )

    run_range = st.sidebar.slider(
        "Run range",
        int(df["run"].min()),
        int(df["run"].max()),
        (int(df["run"].min()), int(df["run"].max()))
    )

    filtered = df[(df["run"] >= run_range[0]) & (df["run"] <= run_range[1])]

    # --- KPIs row ---
    st.subheader("Summary KPIs")
    col1, col2, col3, col4 = st.columns(4)

    col1.metric("Best Objective (max)", f"{filtered['best_obj'].max():.2f}")
    col2.metric("Avg Aircraft Utilization", f"{filtered['aircraft_utilization'].mean():.3f}")
    col3.metric("Avg Passenger Served %", f"{filtered['passenger_demand_served_pct'].mean():.2f}")
    col4.metric("Avg Idle Time", f"{filtered['aircraft_idle_time'].mean():.1f}")

    # --- Line charts ---
    st.subheader("Trends over Runs")

    for metric in selected_metrics:
        fig = px.line(
            filtered,
            x="run",
            y=metric,
            markers=True,
            title=f"{metric} over runs"
        )
        st.plotly_chart(fig, use_container_width=True)

    # --- Correlation heatmap ---
    st.subheader("Metric Correlation")
    corr = filtered[metric_options].corr()

    fig = px.imshow(
        corr,
        text_auto=True,
        aspect="auto",
        title="Correlation Heatmap"
    )
    st.plotly_chart(fig, use_container_width=True)

