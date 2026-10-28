"""Reader for the .vtu/.pvtu/.pvd files written by estran.

Only numpy: the format is the one in src/output.f90 (UInt64 headers,
little endian, uncompressed, cell data), inline base64 or raw appended,
not general VTK.
"""
import base64
import os
import re
import xml.etree.ElementTree as ET
import numpy as np

_DT = {"Float64": "<f8", "Int32": "<i4", "UInt8": "u1"}


def _shape(name, a):
    return a.reshape(-1, 3) if name in ("points", "connectivity") else a


def read_piece(path):
    raw = open(path, "rb").read()
    out = {}
    if b'<AppendedData encoding="raw">' in raw:
        start = raw.index(b'<AppendedData encoding="raw">')
        data = raw.index(b"_", start) + 1
        for m in re.finditer(r"<DataArray ([^>]*)/>", raw[:start].decode()):
            attr = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1)))
            off = data + int(attr["offset"])
            n = int(np.frombuffer(raw, "<u8", 1, off)[0])
            name = attr.get("Name", "points")
            out[name] = _shape(name, np.frombuffer(raw, _DT[attr["type"]], n // np.dtype(_DT[attr["type"]]).itemsize,
                                                   off + 8))
        return out
    for m in re.finditer(rb"<DataArray ([^>]*)>\s*([A-Za-z0-9+/=]*)\s*</DataArray>", raw):
        attr = dict(re.findall(r'(\w+)="([^"]*)"', m.group(1).decode()))
        b = base64.b64decode(m.group(2))
        n = int(np.frombuffer(b, "<u8", 1)[0])
        name = attr.get("Name", "points")
        out[name] = _shape(name, np.frombuffer(b[8:8 + n], _DT[attr["type"]]))
    return out


def read(path):
    """A .pvtu (all pieces merged) or a single .vtu."""
    if path.endswith(".vtu"):
        return read_piece(path)
    root = ET.parse(path).getroot()
    pieces = [read_piece(os.path.join(os.path.dirname(path), p.get("Source")))
              for p in root.iter("Piece")]
    out, shift = {}, 0
    for p in pieces:
        p["connectivity"] = p["connectivity"] + shift
        shift += len(p["points"])
    for k in pieces[0]:
        if k in ("offsets", "types"):
            continue
        out[k] = np.concatenate([p[k] for p in pieces])
    return out


def series(pvd):
    """[(time, path to .pvtu)] from a .pvd collection."""
    root = ET.parse(pvd).getroot()
    d = os.path.dirname(pvd)
    return [(float(e.get("timestep")), os.path.join(d, e.get("file"))) for e in root.iter("DataSet")]


def centroids(f):
    return f["points"][f["connectivity"]][:, :, :2].mean(axis=1)


def areas(f):
    p = f["points"][f["connectivity"]][:, :, :2]
    a, b = p[:, 1] - p[:, 0], p[:, 2] - p[:, 0]
    return 0.5 * (a[:, 0] * b[:, 1] - a[:, 1] * b[:, 0])
