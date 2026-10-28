"""Lake at rest and LeVeque's pulse over the bump.

    python scripts/plot_bump.py        (after estran cases/lake.nml and cases/bump.nml)
"""
import os

import matplotlib.patheffects as pe
import matplotlib.tri as mtri
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

import vtu
from common import BOLD, RES, plt, save

f = vtu.read(os.path.join(RES, "lake", "lake_00002.pvtu"))
print(f"lake at rest after 1 s: max |u| = {np.abs(f['u']).max():.1e} m/s, "
      f"max |v| = {np.abs(f['v']).max():.1e} m/s, max |eta - 1| = {np.abs(f['eta'] - 1).max():.1e} m")

frames = vtu.series(os.path.join(RES, "bump", "bump.pvd"))[1:6]
halo = [pe.withStroke(linewidth=3, foreground="white")]
fig, axes = plt.subplots(3, 2, figsize=(9.2, 7.4), sharex=True, sharey=True, layout="constrained")
f = vtu.read(frames[0][1])
tri = mtri.Triangulation(f["points"][:, 0], f["points"][:, 1], f["connectivity"])
zb = 0.8 * np.exp(-5 * (tri.x - 0.9) ** 2 - 50 * (tri.y - 0.5) ** 2)
# first panel: the bed itself, so the reader knows what the pulse is crossing
ax = axes.flat[0]
sand = LinearSegmentedColormap.from_list("sand", ["#f8f5ee", "#d9c497", "#8a6d3b"])
ax.tricontourf(tri, zb, levels=np.linspace(0, 0.8, 9), cmap=sand)
ax.text(0.03, 0.92, "bed elevation, 0 to 0.8 m", transform=ax.transAxes, va="top", fontsize=9, path_effects=halo)
for ax, (t, path) in zip(axes.flat[1:], frames):
    f = vtu.read(path)
    im = ax.tripcolor(tri, facecolors=1e3 * (f["eta"] - 1), cmap="RdBu_r", vmin=-8, vmax=8, rasterized=True)
    ax.text(0.03, 0.92, f"t = {t:.2f} s", transform=ax.transAxes, va="top", fontsize=9, fontweight=BOLD,
            path_effects=halo)
for ax in axes.flat:
    ax.tricontour(tri, zb, levels=[0.2, 0.4, 0.6], colors="#111827", linewidths=0.5, alpha=0.45)
    ax.set_aspect("equal")
    ax.grid(False)
for ax in axes[:, 0]:
    ax.set_ylabel("y (m)")
for ax in axes[-1]:
    ax.set_xlabel("x (m)")
cb = fig.colorbar(im, ax=axes, orientation="horizontal", shrink=0.45, aspect=40, pad=0.02)
cb.set_label("free surface above 1 m (mm)")
cb.outline.set_visible(False)
fig.suptitle("A 1 cm pulse crossing the bump (LeVeque 1998)", x=0.01, ha="left", fontsize=12, fontweight=BOLD)
save(fig, "bump.png")