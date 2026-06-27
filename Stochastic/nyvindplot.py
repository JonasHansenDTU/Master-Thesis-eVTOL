import numpy as np
import matplotlib.pyplot as plt

# =====================================================
# SETTINGS & FONT STYLING 
# =====================================================
plt.rcParams.update({
    "font.family": "serif",
    "font.size": 13,             
    "axes.titlesize": 26,        
    "axes.labelsize": 16,
    "xtick.labelsize": 16,       
    "ytick.labelsize": 15        
})

# De 12 DMI-sektorer i grader (omregnet til radianer)
compass_degrees = np.array([0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330])
theta = np.radians(compass_degrees)

# DINE ENDELIGE DATA (Opdelt i DMI's 3 kategorier)
light_vinde    = np.array([3.1, 3.1, 3.1, 4.4, 5.8, 5.1, 5.5, 6.7, 7.5, 7.6, 4.9, 2.6]) # 1-9 kts
moderate_vinde = np.array([1.0, 1.0, 1.5, 2.7, 3.5, 2.1, 2.1, 3.7, 5.5, 6.7, 4.1, 0.9]) # 10-19 kts
starke_vinde   = np.array([0.1, 0.0, 0.1, 0.1, 0.3, 0.1, 0.1, 0.3, 0.7, 0.7, 0.5, 0.0]) # >=20 kts

# Bredde på søjlerne
width = np.radians(24) 

# Kvadratisk layout til polarplot
fig, ax = plt.subplots(figsize=(8, 8), subplot_kw=dict(polar=True))

# Sæt Nord (0°) til øverst og tæl MED uret
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

# =====================================================
# PLOTTING (Matchende Folium farvepalette)
# =====================================================
# 1. Lag: Let vind - Lys baggrunds rød fra kortet
bars1 = ax.bar(theta, light_vinde, width=width, color="#F9EAEA", 
               label="Light Wind (1-9 kts)", edgecolor="#641818", linewidth=0.5, alpha=0.9, zorder=2)

# 2. Lag: Moderate vinde - Mellem-rød fra kortet
bars2 = ax.bar(theta, moderate_vinde, width=width, bottom=light_vinde, color="#C95F5F", 
               label="Moderate Wind (10-19 kts)", edgecolor="#641818", linewidth=0.5, alpha=0.9, zorder=3)

# 3. Lag: Stærke vinde - Mørk bordeaux/grænsefarve fra kortet
bars3 = ax.bar(theta, starke_vinde, width=width, bottom=light_vinde + moderate_vinde, color="#641818", 
               label="Strong Wind (>=20 kts)", edgecolor="#641818", linewidth=0.5, alpha=0.95, zorder=4)

# Sæt de 12 kompaslabels på de rigtige pladser
ax.set_xticks(theta)
ax.set_xticklabels(['N', '30°', '60°', 'E', '120°', '150°', 'S', '210°', '240°', 'W', '300°', '330°'])
ax.tick_params(axis='x', pad=15) 

# Justering af procent-cirklerne (vises i en 45 graders vinkel)
ax.set_rlabel_position(45) 
ax.tick_params(axis='y', colors="#641818")

# Tilføj et let baggrundsnet
plt.grid(True, alpha=0.15)

# Titel og legende
plt.title("Operational Wind Rose for\neVTOL Network Scenarios", fontweight="bold", pad=40, color="#333333")
plt.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1, frameon=False, fontsize=18)

# Skab luft i top og bund, så intet bliver beskåret
plt.subplots_adjust(bottom=0.20, top=0.80)

# Gemmer billedet i høj opløsning
plt.savefig('Folium_Matched_Vindrose_Kts.png', bbox_inches='tight', dpi=300)
plt.show()