"""Thacker's paraboloid: free surface across the bowl and error over three
periods.

    python scripts/plot_thacker.py
"""
import os
import numpy as np
import exact
import vtu
from common import AQUA, BLUE, ORANGE, RES, plt, save

frames = vtu.series(os.path.join(RES, "thacker", "thacker.pvd"))
T = exact.PERIOD
pick = [frames[2], frames[4], frames[-1]]          # T/4, T/2, 3T
xx = np.linspace(0, 4, 801)

fig, ax = plt.subplots(figsize=(8, 4.2), layout="constrained")
zb = exact.thacker_bed(xx, 2.0)
ax.fill_between(xx, zb, -0.12, color="#d6c7a1", lw=0)
for (t, path), col in zip(pick, [AQUA, ORANGE, BLUE]):
    f = vtu.read(path)
    c = vtu.centroids(f)
    m = (np.abs(c[:, 1] - 2.0) < 0.03) & (f["h"] > 2e-3)    # leave out the films left by the receding shore
    o = np.argsort(c[m, 0])
    he = exact.thacker(xx, 2.0, t)
    eta = np.where(he > 0, he + zb, np.nan)
    ax.plot(xx, eta, color="#111827", lw=0.8, ls="--")
    ax.plot(c[m, 0][o], f["eta"][m][o], "o", color=col, ms=3, mew=0, label=f"t = {t / T:.2f} T")
ax.plot([], [], color="#111827", ls="--", lw=0.8, label="exact")
ax.set_xlim(0.6, 3.4)
ax.set_ylim(-0.11, 0.105)
ax.set_xlabel("x (m), section y = 2 m")
ax.set_ylabel("elevation (m)")
ax.legend(ncol=4, loc="upper center")
ax.set_title("Thacker's paraboloid, free surface across the bowl", fontsize=11)
save(fig, "thacker.png")

e = np.loadtxt(os.path.join(RES, "thacker", "errors.csv"), delimiter=",", skiprows=1)
fig, ax = plt.subplots(figsize=(6, 3.6), layout="constrained")
ax.fill_between(e[:, 0] / T, 0, 1e3 * e[:, 1], color=BLUE, alpha=0.08, lw=0)
ax.plot(e[:, 0] / T, 1e3 * e[:, 1], "o-", color=BLUE, ms=4, lw=1.5)
ax.set_xlim(0, e[-1, 0] / T)
ax.set_ylim(0, None)
ax.set_title("Mean depth error over three periods (h0 = 100 mm)", fontsize=11)
ax.set_xlabel("time / period")
ax.set_ylabel("mean |h - h_exact| (mm)")
save(fig, "thacker_error.png")

print(f"after 3 periods: mean |h - h_exact| = {1e3 * e[-1, 1]:.2f} mm over the basin, "
      f"max = {1e3 * e[-1, 2]:.1f} mm (at the moving shoreline), h0 = 100 mm")
