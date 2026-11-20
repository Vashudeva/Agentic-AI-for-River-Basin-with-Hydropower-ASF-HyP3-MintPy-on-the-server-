"""Self-check before (or without) the long processing.

Writes MintPy-format geocoded files for a variant called 'synthetic' that contain the
SARscape LOS histories themselves (each SARscape point fills the grid cell it falls in).
Running the products and compare steps on it must then give:
  * ascending / descending: zero differences (tests sampling, re-referencing, statistics);
  * vertical / east: the agreement between kgproc's decomposition and SARscape's own
    decomposition of the same LOS data (tests unit vectors, signs, date interpolation).

    python tests/synthetic_from_sarscape.py -c config.toml
    python -m kgproc products --variant synthetic -c config.toml
    python -m kgproc compare  --variant synthetic -c config.toml
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kgproc.config import TRACKS, load                       # noqa: E402
from kgproc.products import read_sarscape                     # noqa: E402


def mintpy_az_from_sarscape(alos: np.ndarray) -> np.ndarray:
    """SARscape ALOS (about +80 asc / -79 dsc) -> MintPy LOS azimuth (about +100 / -101)."""
    return np.sign(alos) * (180.0 - np.abs(alos))


def write(path: Path, datasets: dict, attrs: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with h5py.File(path, "w") as f:
        for k, v in datasets.items():
            f.create_dataset(k, data=v)
        for k, v in attrs.items():
            f.attrs[k] = str(v)


def main(argv=None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--variant", default="synthetic")
    a = ap.parse_args(argv)
    cfg = load(a.config)
    S, N, W, E = cfg.mintpy["geocode_snwe"]
    step = float(cfg.mintpy.get("geocode_step_deg", 0.000125))
    rows, cols = int(round((N - S) / step)), int(round((E - W) / step))
    lat_ref, lon_ref = cfg.mintpy["reference_lalo"]
    base = {"Y_FIRST": N, "X_FIRST": W, "Y_STEP": -step, "X_STEP": step, "LENGTH": rows, "WIDTH": cols,
            "REF_LAT": lat_ref, "REF_LON": lon_ref, "PROCESSOR": "synthetic", "UNIT": "m"}
    for key in TRACKS:
        fr, dcols, lon, lat, src = read_sarscape(cfg.sarscape_products[key])
        D = fr[dcols].apply(pd.to_numeric, errors="coerce").to_numpy(np.float32) / 1000.0     # m
        i = np.floor((lat - N) / -step).astype(int)
        j = np.floor((lon - W) / step).astype(int)
        ok = (i >= 0) & (i < rows) & (j >= 0) & (j < cols)
        ts = np.full((len(dcols), rows, cols), np.nan, np.float32)
        ts[:, i[ok], j[ok]] = D[ok].T
        tcoh = np.zeros((rows, cols), np.float32)
        tcoh[i[ok], j[ok]] = 0.9
        geom = {}
        for name, field, conv in (("height", "Z", None), ("incidenceAngle", "ILOS", None),
                                  ("azimuthAngle", "ALOS", mintpy_az_from_sarscape)):
            g = np.full((rows, cols), np.nan, np.float32)
            v = pd.to_numeric(fr[field], errors="coerce").to_numpy(float)
            g[i[ok], j[ok]] = conv(v[ok]) if conv else v[ok]
            fill = np.nanmedian(g)
            geom[name] = np.where(np.isfinite(g), g, fill).astype(np.float32)
        vel = np.full((rows, cols), np.nan, np.float32)
        vel[i[ok], j[ok]] = pd.to_numeric(fr["velocity"], errors="coerce").to_numpy(float)[ok] / 1000.0
        out = cfg.geo_dir(key, a.variant)
        dates = np.array([c[2:] for c in dcols], dtype="S8")
        write(out / "geo_timeseries.h5", {"timeseries": ts, "date": dates,
                                          "bperp": np.zeros(len(dcols), np.float32)},
              {**base, "FILE_TYPE": "timeseries"})
        write(out / "geo_temporalCoherence.h5", {"temporalCoherence": tcoh}, {**base, "FILE_TYPE": "temporalCoherence"})
        write(out / "geo_geometryRadar.h5", geom, {**base, "FILE_TYPE": "geometry"})
        write(out / "geo_velocity.h5", {"velocity": vel}, {**base, "FILE_TYPE": "velocity", "UNIT": "m/year"})
        print(f"{key}: {ok.sum()} SARscape points written to {out} ({rows} x {cols} grid, {len(dcols)} dates)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
