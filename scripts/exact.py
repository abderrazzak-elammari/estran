"""Exact solutions, same formulas as src/cases.f90"""
import numpy as np

g = 9.81

# Thacker (1981), planar surface in a paraboloid (SWASHES notation)
A, H0, ETA, L = 1.0, 0.1, 0.5, 4.0
OMEGA = np.sqrt(2 * g * H0) / A
PERIOD = 2 * np.pi / OMEGA


def thacker_bed(x, y):
    return -H0 * (1 - ((x - L / 2) ** 2 + (y - L / 2) ** 2) / A**2)


def thacker(x, y, t):
    h = ETA * H0 / A**2 * (2 * (x - L / 2) * np.cos(OMEGA * t) + 2 * (y - L / 2) * np.sin(OMEGA * t) - ETA)
    return np.maximum(h - thacker_bed(x, y), 0.0)


def dambreak(x, t, hl=1.0, hr=0.1, x0=25.0):
    """Stoker (hr > 0) or Ritter (hr = 0)"""
    xs = np.asarray(x, float) - x0
    cl = np.sqrt(g * hl)
    if hr <= 0:
        h = np.where(xs < -cl * t, hl, (2 * cl - xs / t) ** 2 / (9 * g))
        return np.where(xs > 2 * cl * t, 0.0, h)
    lo, hi = hr, hl
    for _ in range(200):
        hm = 0.5 * (lo + hi)
        f = 2 * (cl - np.sqrt(g * hm)) - (hm - hr) * np.sqrt(0.5 * g * (hm + hr) / (hm * hr))
        lo, hi = (hm, hi) if f > 0 else (lo, hm)
    cm = np.sqrt(g * hm)
    um = 2 * (cl - cm)
    s = hm * um / (hm - hr)
    return np.select([xs < -cl * t, xs < (um - cm) * t, xs < s * t],
                     [hl, (2 * cl - xs / t) ** 2 / (9 * g), hm], hr)
