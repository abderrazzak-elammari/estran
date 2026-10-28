import os
import subprocess
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from matplotlib import font_manager

FONT = "Inter" if "Inter" in {f.name for f in font_manager.fontManager.ttflist} else "DejaVu Sans"
BOLD = 600 if FONT == "Inter" else "bold"      

plt.rcParams.update({"font.family": FONT,"font.size": 10, "axes.titlesize": 11, "axes.titleweight": BOLD,
    "axes.titlelocation": "left","axes.titlepad": 8,"axes.labelcolor": "#374151","axes.edgecolor": "#9ca3af","axes.linewidth": 0.8,
    "axes.spines.top": False,"axes.spines.right": False,"axes.grid": True,"axes.axisbelow": True,"grid.color": "#e5e7eb","grid.linewidth": 0.6,
    "xtick.color": "#4b5563","ytick.color": "#4b5563","xtick.labelsize": 9,"ytick.labelsize": 9,"legend.frameon": False,"legend.fontsize": 9,"lines.linewidth": 1.8,"figure.facecolor": "white",
})



ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
RES = os.path.join(ROOT, "results")

BLUE, ORANGE, AQUA, RED, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#d03b3b", "#6b7280"


def run(nml, np_=2, overrides=None, quiet=True):
    """Run estran on a namelist, optionally with keys replaced (key=value strings)."""
    if overrides:
        txt = open(os.path.join(ROOT, nml)).read().rstrip().rstrip("/")
        txt += "\n" + "\n".join("  " + o for o in overrides) + "\n/\n"
        nml = os.path.join(RES, "tmp_" + os.path.basename(nml))
        os.makedirs(RES, exist_ok=True)
        open(nml, "w").write(txt)
    cmd = os.environ.get("MPIRUN", "mpirun").split() + ["-np", str(np_), os.path.join(ROOT, "estran"), nml]
    out = subprocess.run(cmd, cwd=ROOT, capture_output=quiet, text=True)
    if out.returncode:
        sys.exit(out.stderr or "estran failed")
    return out.stdout


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, color="#e5e7eb", lw=0.6)
    ax.set_axisbelow(True)


def save(fig, name):
    os.makedirs(FIG, exist_ok=True)
    fig.savefig(os.path.join(FIG, name), dpi=200, bbox_inches="tight", facecolor="white")
    print("figures/" + name)
