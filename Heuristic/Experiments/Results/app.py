import pandas as pd
from Dashborad import run_dashboard

df = pd.read_csv("Heuristic/Experiments/Results/best_solution_kpis.csv")
run_dashboard(df)