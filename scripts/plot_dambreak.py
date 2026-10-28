"""Stoker and Ritter dam breaks against the exact solutions, and the L1
convergence on the four meshes.

    python scripts/plot_dambreak.py
"""
import os
import numpy as np
import exact
import vtu
from common import BLUE, BOLD, GREY, RES, plt, save

def centreline(path, band=0.06):
    f = vtu.read(path)
    c = vtu.centroids(f)
    m = np.abs(c[:, 1] - 1.0) < band
    o = np.argsort(c[m, 0])
    return c[m, 0][o], f["h"][m][o]


x = np.linspace(0, 50, 4001)
ZOOM = {"dambreak_4": ((33.4, 35.4), (0.05, 0.45), [0.64, 0.42, 0.33, 0.45]),
        "dambreak_dry": ((39.0, 45.0), (-0.002, 0.03), [0.6, 0.42, 0.37, 0.45])}
fig, axes = plt.subplots(1, 2, figsize=(12, 4.2), sharey=True, layout="constrained")
for ax, (d, hr, title) in zip(axes, [("dambreak_4", 0.1, "Stoker, wet bed (h_r = 0.1 m)"),
                                     ("dambreak_dry", 0.0, "Ritter, dry bed")]):
    pvd = [p for p in os.listdir(os.path.join(RES, d)) if p.endswith(".pvd")][0]
    t, path = vtu.series(os.path.join(RES, d, pvd))[-1]
    xs, h = centreline(path)
    he = exact.dambreak(x, t, hr=hr)
    zx, zy, box = ZOOM[d]
    ins = ax.inset_axes(box)
    for a in (ax, ins):
        a.fill_between(x, 0, he, color="#dbeafe", lw=0)
        a.plot(x, he, color="#111827", lw=0.9, label="exact")
        a.plot(xs, h, "o", color=BLUE, ms=2.2 if a is ax else 3.5, mew=0, label="estran, cells on the centreline")
    ins.set_xlim(*zx)
    ins.set_ylim(*zy)
    ins.tick_params(labelsize=7)
    ins.grid(False)
    for sp in ins.spines.values():
        sp.set_visible(True)
        sp.set_color("#9ca3af")
    ax.indicate_inset_zoom(ins, edgecolor="#6b7280", alpha=0.8)
    ax.set_title(f"{title}, t = {t:.0f} s", fontsize=10)
    ax.set_xlabel("x (m)")
    ax.set_xlim(10, 50)
    ax.set_ylim(-0.03, 1.08)
axes[0].set_ylabel("water depth h (m)")
axes[1].legend(loc="lower left")
fig.suptitle("Dam break on a flat bed, estran against the exact solution", x=0.01, ha="left", fontsize=12,
             fontweight=BOLD)
save(fig, "dambreak.png")

lc = np.array([0.5, 0.25, 0.125, 0.0625])
err = np.array([np.loadtxt(os.path.join(RES, f"dambreak_{k}", "errors.csv"), delimiter=",", skiprows=1)[-1, 1]
                for k in range(1, 5)])
cells = [len(vtu.read(vtu.series(os.path.join(RES, f"dambreak_{k}", "dambreak.pvd"))[-1][1])["h"]) for k in range(1, 5)]
print("mesh size (m)   cells    L1 error (m)   order")
for k in range(4):
    order = "" if k == 0 else f"{np.log(err[k-1] / err[k]) / np.log(2):.2f}"
    print(f"{lc[k]:12.4f}  {cells[k]:7d}  {err[k]:13.3e}   {order}")

fig, ax = plt.subplots(figsize=(5.4, 4), layout="constrained")
ax.loglog(lc, err, "o-", color=BLUE, lw=2, ms=6, label="estran")
ax.loglog(lc, err[0] * lc / lc[0], "--", color=GREY, lw=1, label="first order")
for k in range(1, 4):
    p = np.log(err[k - 1] / err[k]) / np.log(2)
    ax.annotate(f"{p:.2f}", (np.sqrt(lc[k - 1] * lc[k]), np.sqrt(err[k - 1] * err[k])), xytext=(6, -12),
                textcoords="offset points", fontsize=8, color=BLUE)
ax.set_xticks(lc, [f"{v:g}" for v in lc])
ax.set_yticks([3e-4, 5e-4, 1e-3, 2e-3], ["3e-4", "5e-4", "1e-3", "2e-3"])
ax.minorticks_off()
ax.set_xlabel("mesh size (m)")
ax.set_ylabel("L1 error on h at t = 3 s (m)")
ax.legend(loc="upper left")
ax.set_title("L1 convergence, Stoker case (observed order in blue)", fontsize=10)
save(fig, "dambreak_convergence.png")
