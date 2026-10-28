"""Strong scaling on the island mesh: 2 s of case A, no output.

    python scripts/scaling.py 1 2 4 8
"""
import os
import sys
import numpy as np
from common import BLUE, BOLD, GREY, ORANGE, RES, plt, run, save

replot = "--plot" in sys.argv
ranks = [int(a) for a in sys.argv[1:] if a != "--plot"] or [1, 2, 4]
rows = []
for n in ranks:
    d = f"scaling_{n}"
    if not replot:
        run("cases/island_a.nml", n, [f"out_dir = 'results/{d}'", "t_end = 2.0", "dt_out = 0", "ng = 0"])
    r = np.loadtxt(os.path.join(RES, d, "timing.csv"), delimiter=",", skiprows=1)
    rows.append(r)
rows = np.array(rows)
n, wall, comm = rows[:, 0], rows[:, 3], rows[:, 4]
speedup = wall[0] / wall
print("ranks   wall (s)   speed-up   efficiency   comm on rank 0")
for k in range(len(n)):
    print(f"{int(n[k]):5d}  {wall[k]:9.2f}  {speedup[k]:9.2f}  {speedup[k] / n[k] * n[0]:10.0%}  {comm[k] / wall[k]:14.0%}")

eff = speedup / n * n[0]
fig, (ax, bx) = plt.subplots(1, 2, figsize=(11, 4.2), layout="constrained")
ax.plot(n, n / n[0], "--", color=GREY, lw=1, label="ideal")
ax.plot(n, speedup, "o-", color=BLUE, lw=2, ms=7, label="estran")
for k in range(1, len(n)):
    ax.annotate(f"{eff[k]:.0%}", (n[k], speedup[k]), xytext=(0, -16), textcoords="offset points", ha="center",
                fontsize=8, color=BLUE)
ax.set_xlabel("MPI ranks")
ax.set_ylabel("speed-up")
ax.set_xticks(n)
ax.set_ylim(0, 1.08 * n[-1] / n[0])
ax.legend(loc="upper left")
ax.set_title("Speed-up, parallel efficiency under each point", fontsize=10)

x = np.arange(len(n))
bx.bar(x, wall - comm, 0.6, color=BLUE, label="computation")
bx.bar(x, comm, 0.6, bottom=wall - comm, color=ORANGE, label="communication and waiting, rank 0")
for k in range(len(n)):
    bx.text(x[k], wall[k], f"{wall[k]:.1f} s", ha="center", va="bottom", fontsize=8, color="#374151")
bx.set_xticks(x, [str(int(v)) for v in n])
bx.set_ylim(0, 1.12 * wall.max())
bx.set_xlabel("MPI ranks")
bx.set_ylabel("wall time (s)")
bx.grid(axis="x", visible=False)
bx.legend(loc="upper right")
bx.set_title("Where the time goes", fontsize=10)
fig.suptitle(f"Strong scaling, conical island, {int(rows[0, 1])} cells, {int(rows[0, 2])} time steps", x=0.01,
             ha="left", fontsize=12, fontweight=BOLD)
save(fig, "scaling.png")
