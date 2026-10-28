"""How the island mesh is cut between ranks.

    python scripts/plot_partition.py [nranks]       (default 8, oversubscribes if needed)
"""
import os
import sys
import matplotlib.patheffects as pe
import matplotlib.tri as mtri
import numpy as np
from matplotlib.colors import ListedColormap
import vtu
from common import RES, plt, run, save

NP = int(sys.argv[1]) if len(sys.argv) > 1 else 8
COLORS = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]

d = f"partition_{NP}"
run("cases/island_a.nml", NP, [f"out_dir = 'results/{d}'", "t_end = 1e-3", "dt_out = 1e-3", "ng = 0"])
f = vtu.read(vtu.series(os.path.join(RES, d, "island.pvd"))[0][1])
tri = mtri.Triangulation(f["points"][:, 0], f["points"][:, 1], f["connectivity"])
c = vtu.centroids(f)

fig, ax = plt.subplots(figsize=(7.5, 8))
ax.tripcolor(tri, facecolors=f["rank"].astype(float), cmap=ListedColormap(COLORS[:NP] if NP <= 8 else None),
             vmin=-0.5, vmax=NP - 0.5, edgecolors="white", linewidth=0.05, rasterized=True)
for r in range(NP):
    m = f["rank"] == r
    ax.text(*np.median(c[m], axis=0), f"rank {r}\n{m.sum()} cells", ha="center", va="center", fontsize=8,
            color="#111827", path_effects=[pe.withStroke(linewidth=2.5, foreground="white")])
th = np.linspace(0, 2 * np.pi, 200)
for rad in (3.6, 1.1):
    ax.plot(12.96 + rad * np.cos(th), 13.80 + rad * np.sin(th), color="white", lw=0.8)
ax.set_aspect("equal")
ax.set_xlabel("x (m)")
ax.set_ylabel("y (m)")
ax.set_title(f"{len(f['h'])} triangles on {NP} ranks: equal counts, so small parts where the mesh is fine",
             fontsize=10)
save(fig, f"partition_{NP}.png")
