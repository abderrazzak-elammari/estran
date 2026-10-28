"""Conical island against the Briggs et al. (1995) measurements: surface
elevation at the gauges and maximum run-up around the island.

    python scripts/plot_island.py
"""
import os
import matplotlib.patheffects as pe
import numpy as np
import vtu
from common import BLUE, BOLD, GREY, RES, ROOT, plt, save

DATA = os.path.join(ROOT, "data", "briggs1995")
XC, YC = 12.96, 13.80
NAMES = ["1", "2", "3", "4", "6", "9", "16", "22"]
CASES = [("a", "ts2a.txt", "run2a.txt", "A, H/d = 0.045"), ("c", "ts2cnew1.txt", "run2c.txt", "C, H/d = 0.181")]


def gauges(case, ts):
    sim = np.loadtxt(os.path.join(RES, f"island_{case}", "gauges.csv"), delimiter=",", skiprows=1)
    exp = np.loadtxt(os.path.join(DATA, ts), skiprows=7)
    # the gauges do not read zero in still water (up to 3 mm off): take the
    # mean of the first two seconds, before the wave arrives, as zero
    exp[:, 1:] -= exp[exp[:, 0] < 22, 1:].mean(axis=0)
    # the wavemaker clock and ours differ by a constant: line up the crests
    # at gauge 2, the middle of the incident line
    shift = exp[np.argmax(exp[:, 2]), 0] - sim[np.argmax(sim[:, 2]), 0]
    return sim, exp, shift


def runup(case, dtheta=5.0, h_wet=1e-3):
    f = vtu.read(vtu.series(os.path.join(RES, f"island_{case}", "island.pvd"))[-1][1])
    c = vtu.centroids(f)
    dx, dy = c[:, 0] - XC, c[:, 1] - YC
    th = np.degrees(np.arctan2(dx, -dy)) % 360       # 0 deg towards -y, 90 deg on the lee side
    m = (np.hypot(dx, dy) < 3.6) & (f["hmax"] > h_wet)
    bins = np.arange(0, 360 + dtheta, dtheta)
    k = np.digitize(th[m], bins) - 1
    z = np.full(len(bins) - 1, np.nan)
    for i in range(len(z)):
        if np.any(k == i):
            z[i] = f["b"][m][k == i].max()
    return 0.5 * (bins[1:] + bins[:-1]), z


ROLE = ["incident", "incident", "incident", "incident", "toe, front", "slope, front", "slope, side", "slope, lee"]
GX = {"a": 5.76, "c": 7.56}          # the incident line was moved between the two cases
GPOS = [(None, 16.05), (None, 14.55), (None, 13.05), (None, 11.55), (9.36, 13.80), (10.36, 13.80), (12.96, 11.22),
        (15.56, 13.80)]
# label offset (m) and alignment for each gauge on the map
GLAB = [(-0.35, 0, "right")] * 4 + [(0, 0.5, "center"), (0, 0.5, "center"), (0, -0.55, "center"), (0.35, 0, "left")]
WATER, SAND, SAND_DARK, INK = "#eef5fc", "#e9dcb9", "#c2a46a", "#111827"


def gauge_map(ax, case):
    """Plan view around the island: toe, shoreline at rest, crest and gauges."""
    th = np.linspace(0, 2 * np.pi, 200)
    ax.set_facecolor(WATER)
    for r, kw in [(3.6, dict(color=SAND)), (2.32, dict(color=SAND_DARK)), (1.1, dict(color="#a88a4f"))]:
        ax.fill(XC + r * np.cos(th), YC + r * np.sin(th), lw=0, **kw)
    for j, (x, y) in enumerate(GPOS):
        x = GX[case] if x is None else x
        ax.plot(x, y, "o", ms=6, mfc=BLUE, mec="white", mew=1)
        dx, dy, ha = GLAB[j]
        ax.text(x + dx, y + dy, NAMES[j], ha=ha, va="center", fontsize=8, color=INK)
    ax.annotate("", xy=(7.3, 10.1), xytext=(4.9, 10.1), arrowprops=dict(arrowstyle="-|>", color="#374151", lw=1.3))
    ax.text(4.9, 10.55, "incident wave", fontsize=8, color="#374151")
    ax.text(XC, YC + 3.0, "shoreline at rest", ha="center", fontsize=7, color="#6b5428")
    ax.set_xlim(4.5, 17.5)
    ax.set_ylim(9.5, 18.2)
    ax.set_aspect("equal")
    ax.grid(False)
    for sp in ax.spines.values():
        sp.set_visible(True)
    ax.tick_params(labelsize=8)
    ax.set_xlabel("x (m)", fontsize=9)
    ax.set_ylabel("y (m)", fontsize=9)
    ax.set_title("Gauge layout", fontsize=10)


def peak_panel(ax, exp, sim):
    """Highest elevation at each gauge, measured against computed."""
    y = np.arange(8)[::-1]
    pm, ps = 100 * exp[:, 1:9].max(axis=0), 100 * sim[:, 1:9].max(axis=0)
    ax.hlines(y, np.minimum(pm, ps), np.maximum(pm, ps), color="#d1d5db", lw=2.5)
    ax.plot(pm, y, "o", ms=6, mfc="white", mec=GREY, mew=1.4, label="measured")
    ax.plot(ps, y, "o", ms=6, color=BLUE, label="estran")
    ax.set_yticks(y, [f"gauge {n}" for n in NAMES], fontsize=8)
    ax.grid(axis="y", visible=False)
    ax.set_xlim(0, 1.15 * max(pm.max(), ps.max()))
    ax.set_xlabel("peak elevation (cm)", fontsize=9)
    ax.legend(loc="upper left", fontsize=8, handletextpad=0.3)
    ax.set_title("Highest crest at each gauge", fontsize=10)


for case, ts, run, label in CASES:
    sim, exp, shift = gauges(case, ts)
    fig = plt.figure(figsize=(14, 8.5), layout="constrained")
    gs = fig.add_gridspec(4, 3, width_ratios=[1, 1, 0.72])
    axes = np.array([[fig.add_subplot(gs[i, j]) for j in range(2)] for i in range(4)])
    for ax in axes.flat[1:]:
        ax.sharex(axes[0, 0])
    for j, ax in enumerate(axes.T.flat):
        ax.plot(exp[:, 0], 100 * exp[:, j + 1], color=GREY, lw=1.1, label="measured")
        ax.plot(sim[:, 0] + shift, 100 * sim[:, j + 1], color=BLUE, lw=1.6, label="estran")
        ax.text(0.015, 0.9, f"gauge {NAMES[j]}", transform=ax.transAxes, fontsize=9, fontweight=BOLD, va="top")
        ax.text(0.015, 0.74, ROLE[j], transform=ax.transAxes, fontsize=8, color=GREY, va="top")
        ax.set_xlim(sim[0, 0] + shift + 2, sim[-1, 0] + shift)
        if ax not in axes[-1]:
            ax.tick_params(labelbottom=False)
    for ax in axes[-1]:
        ax.set_xlabel("time (s), experiment clock")
    axes[0, 1].legend(loc="upper right", ncol=2)
    fig.supylabel("surface elevation (cm)", fontsize=10, color="#374151")
    gauge_map(fig.add_subplot(gs[:2, 2]), case)
    peak_panel(fig.add_subplot(gs[2:, 2]), exp, sim)
    fig.suptitle(f"Conical island, case {label}: surface elevation at the eight gauges", x=0.01, ha="left",
                 fontsize=12, fontweight=BOLD)
    save(fig, f"island_gauges_{case}.png")

    peaks = [(NAMES[j], 100 * exp[:, j + 1].max(), 100 * sim[:, j + 1].max()) for j in range(8)]
    print(f"case {case.upper()}, clock shift {shift:.2f} s. max eta (cm), measured / computed:")
    print("   " + "   ".join(f"g{n}: {a:.2f}/{b:.2f}" for n, a, b in peaks))

# run-up drawn around the island as seen from above: 0 deg points to -y, the
# wave comes from the left (270 deg), the lee side is on the right (90 deg)
halo = [pe.withStroke(linewidth=3, foreground="white")]
fig, axes = plt.subplots(1, 2, figsize=(12, 5.8), subplot_kw=dict(projection="polar"), layout="constrained")
for ax, (case, ts, run, label) in zip(axes, CASES):
    meas = np.loadtxt(os.path.join(DATA, run), skiprows=10)
    th, z = runup(case)
    t = np.radians(np.append(th, th[0] + 360))
    r = 100 * np.append(z, z[0])
    ax.set_theta_zero_location("S")
    ax.set_theta_direction(1)
    ax.fill(t, r, color=BLUE, alpha=0.10, lw=0)
    ax.plot(t, r, color=BLUE, lw=2, label="estran")
    ax.plot(np.radians(meas[:, 1]), meas[:, 2], "o", ms=5.5, mfc="white", mec=INK, mew=1.2, label="measured")
    ax.set_thetagrids([0, 90, 180, 270], ["0°", "90°\nlee side", "180°", "270°\nfacing the wave"], fontsize=8,
                      color="#4b5563")
    ax.tick_params(axis="x", pad=10)
    rmax = 1.15 * max(np.nanmax(r), meas[:, 2].max())
    ticks = np.arange(1, 6) * (2 if case == "c" else 1) * (2.5 if case == "c" else 1)
    ticks = ticks[ticks < rmax]
    ax.set_ylim(0, rmax)
    ax.set_yticks(ticks, [f"{v:g}" for v in ticks[:-1]] + [f"{ticks[-1]:g} cm"], fontsize=8, color="#6b7280")
    for lab in ax.get_yticklabels():
        lab.set_path_effects(halo)
    ax.set_rlabel_position(150)
    ax.grid(color="#e5e7eb", lw=0.7)
    ax.spines["polar"].set_color("#9ca3af")
    ax.set_title(f"Case {label}", pad=16)
    err = np.interp(meas[:, 1], th, 100 * z, period=360) - meas[:, 2]
    print(f"case {case.upper()} run-up: mean |error| {np.mean(np.abs(err)):.2f} cm on "
          f"{len(meas)} points, measured range {meas[:, 2].min():.1f}-{meas[:, 2].max():.1f} cm")
h, l = axes[0].get_legend_handles_labels()
fig.legend(h, l, loc="upper right", ncol=2)
fig.suptitle("Maximum run-up around the island", x=0.01, ha="left", fontsize=12, fontweight=BOLD)
save(fig, "island_runup.png")

# sensitivity of the case A run-up, when those runs exist (see README)
meas = np.loadtxt(os.path.join(DATA, "run2a.txt"), skiprows=10)
for d, what in [("a", "7 cm cells on the island"), ("a_fine", "3.5 cm cells"), ("a_manning", "7 cm cells, Manning n = 0.012")]:
    if not os.path.exists(os.path.join(RES, f"island_{d}", "island.pvd")):
        continue
    th, z = runup(d)
    err = np.interp(meas[:, 1], th, 100 * z, period=360) - meas[:, 2]
    print(f"case A, {what:32s} run-up error: mean {err.mean():+.2f} cm, mean |error| {np.abs(err).mean():.2f} cm")
