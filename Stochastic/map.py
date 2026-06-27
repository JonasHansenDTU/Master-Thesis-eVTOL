import folium

m = folium.Map(location=[55.6, 10.5], zoom_start=8, tiles="CartoDB positron")

locations = [
    ("Billund",    55.740278, 9.151667,  "darkred"),
    ("Kalundborg", 55.7003,  11.2500,    "red"),
    ("Odense",     55.4767,  10.3308,    "lightred"),
]

for name, lat, lon, color in locations:
    folium.Marker(
        location=[lat, lon],
        popup=name,
        tooltip=name,
        icon=folium.Icon(color=color, icon="circle"),
    ).add_to(m)

m.save("/Users/asbb/Desktop/Speciale/Master-Thesis-eVTOL/Stochastic/map.html")
print("Saved!")