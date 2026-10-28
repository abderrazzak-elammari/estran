"""Regression checks, run by `make test`. Short runs, a minute in total.

NP ranks for the parallel ones (default 4), MPIRUN to change the launcher.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import exact  # noqa: E402
import vtu  # noqa: E402
from common import RES, run  # noqa: E402

NP = int(os.environ.get("NP", 4))
failed = []


def check(name, ok, detail):
    print(f"{'ok  ' if ok else 'FAIL'} {name:34s} {detail}")
    if not ok:
        failed.append(name)


def last(d, prefix):
    return vtu.read(vtu.series(os.path.join(RES, d, prefix + ".pvd"))[-1][1])


def sorted_fields(f):
    c = vtu.centroids(f)
    o = np.lexsort((c[:, 1], c[:, 0]))
    return np.stack([f["h"][o], f["u"][o], f["v"][o]])


run("cases/lake.nml", NP, ["out_dir = 'results/test_lake'", "t_end = 0.2", "dt_out = 0.2"])
f = last("test_lake", "lake")
du = max(np.abs(f["u"]).max(), np.abs(f["v"]).max())
deta = np.abs(f["eta"] - 1).max()
check("lake at rest stays at rest", du < 1e-10 and deta < 1e-10, f"max|u| {du:.1e} m/s, max|eta-1| {deta:.1e} m")

a = []
for n in (1, NP):
    run("cases/bump.nml", n, [f"out_dir = 'results/test_bump_{n}'", "t_end = 0.1", "dt_out = 0.1"])
    a.append(sorted_fields(last(f"test_bump_{n}", "bump")))
check(f"1 and {NP} ranks give the same bits", np.array_equal(a[0], a[1]),
      f"max difference {np.abs(a[0] - a[1]).max():.1e}")

run("cases/thacker.nml", NP, ["out_dir = 'results/test_thacker'", f"t_end = {exact.PERIOD / 2}",
                               f"dt_out = {exact.PERIOD / 2}"])
v = np.loadtxt(os.path.join(RES, "test_thacker", "diag.csv"), delimiter=",", skiprows=1)[:, 3]
e = np.loadtxt(os.path.join(RES, "test_thacker", "errors.csv"), delimiter=",", skiprows=1)[-1]
check("closed basin keeps its volume", abs(v[-1] / v[0] - 1) < 1e-12, f"dV/V {v[-1] / v[0] - 1:.1e}")
check("Thacker after half a period", e[1] < 2e-4, f"mean |h - h_exact| {e[1]:.2e} m")

run("cases/dambreak.nml", NP, ["mesh_file = 'meshes/dambreak_2.msh'", "out_dir = 'results/test_dambreak'"])
e = np.loadtxt(os.path.join(RES, "test_dambreak", "errors.csv"), delimiter=",", skiprows=1)[-1]
check("Stoker dam break, 4000 cells", e[1] < 1.2e-3, f"mean |h - h_exact| {e[1]:.2e} m")

run("cases/dambreak_dry.nml", NP, ["mesh_file = 'meshes/dambreak_2.msh'", "out_dir = 'results/test_ritter'"])
f = last("test_ritter", "ritter")
check("dry bed stays non-negative", f["h"].min() >= 0, f"min h {f['h'].min():.1e} m")

print(f"\n{len(failed)} failed" if failed else "\nall passed")
sys.exit(1 if failed else 0)
