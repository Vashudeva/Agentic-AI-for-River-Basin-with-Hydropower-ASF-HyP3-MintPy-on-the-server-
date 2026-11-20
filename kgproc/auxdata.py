"""Step 'aux': precise orbits (sentineleof) and DEM (sardem, Copernicus GLO-30).

ISCE2 needs a DEM in WGS84 ellipsoidal heights covering every burst that touches the
window, i.e. the full width of the sub-swath. sardem converts the Copernicus heights
(EGM2008) to the ellipsoid and writes the ISCE2 .xml.
"""
from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path

import numpy as np

from .config import Config
from .inventory import read_safe
from .util import StepError, run, tool_env

EOF_NAME = re.compile(r"(S1[ABCD])_OPER_AUX_(POEORB|RESORB)_OPOD_\d{8}T\d{6}_V(\d{8}T\d{6})_(\d{8}T\d{6})\.EOF$")


def _inventory(cfg: Config, key: str) -> dict:
    f = cfg.track_dir(key) / "inventory" / "inventory.json"
    if not f.is_file():
        raise StepError(f"Run the inventory step first ({f} missing)")
    return json.loads(f.read_text())


def orbits(cfg: Config, key: str) -> Path:
    inv = _inventory(cfg, key)
    slc = cfg.track_dir(key) / "slc"
    out = cfg.track_dir(key) / "orbits"
    out.mkdir(parents=True, exist_ok=True)
    log = cfg.track_dir(key) / "logs" / "aux.log"
    cmd = ["eof", "--search-path", str(slc), "--save-dir", str(out)]
    if cfg.isce.get("orbit_source", "asf") == "asf":
        cmd.append("--force-asf")
    run(cmd, log)

    # every acquisition must be inside the validity window of an orbit file
    files = []
    for f in out.glob("*.EOF"):
        m = EOF_NAME.search(f.name)
        if m:
            files.append((m.group(1), m.group(2), datetime.strptime(m.group(3), "%Y%m%dT%H%M%S"),
                          datetime.strptime(m.group(4), "%Y%m%dT%H%M%S")))
    lacking, restituted = [], []
    for row in inv["dates"]:
        for fname in row["files"]:
            name = Path(fname).name
            mission, start = name[:3], datetime.strptime(name[17:32], "%Y%m%dT%H%M%S")
            ok = [kind for (ms, kind, a, b) in files if ms == mission and a <= start <= b]
            if not ok:
                lacking.append(name)
            elif "POEORB" not in ok:
                restituted.append(name)
    if lacking:
        raise StepError(f"No orbit file for {len(lacking)} scene(s), e.g. {lacking[0]}. See {log}")
    if restituted:
        print(f"  note: {len(restituted)} scene(s) only have restituted orbits")
    print(f"  orbits ok: {len(files)} file(s) in {out}")
    return out


def dem(cfg: Config) -> Path:
    """One DEM for both tracks, covering the footprints of all scenes near the window."""
    user = str(cfg.isce.get("dem_file") or "").strip()
    out_dir = cfg.work_dir / "dem"
    out_dir.mkdir(parents=True, exist_ok=True)
    log = cfg.work_dir / "logs" / "dem.log"
    if user:
        return _user_dem(Path(user), out_dir, log)

    target = out_dir / "elevation.dem"
    if target.is_file() and target.with_suffix(".dem.xml").is_file():
        _check_dem(target)
        return target
    lons, S, N = [], 90.0, -90.0
    margin = float(cfg.isce.get("dem_margin_deg", 0.3))
    for key in cfg.tracks:
        inv = _inventory(cfg, key)
        s, n, w, e = cfg.tracks[key].window_snwe
        S, N = min(S, s - margin), max(N, n + margin)
        lons += [w - margin, e + margin]
        # the scene footprint gives the full swath width; it repeats along a track,
        # so the first readable scene of each track is enough
        first = next((f for row in inv["dates"] for f in row["files"]), None)
        if first:
            fp = read_safe(Path(first))["footprint"]
            if fp:
                lons += [p[0] for p in fp]
    W, E = min(lons) - 0.05, max(lons) + 0.05
    print(f"  DEM box: W {W:.3f} S {S:.3f} E {E:.3f} N {N:.3f}")
    run(["sardem", "--bbox", f"{W:.4f}", f"{S:.4f}", f"{E:.4f}", f"{N:.4f}", "--data-source", "COP",
         "-isce", "--output-format", "ENVI", "-o", str(target)], log, cwd=out_dir)
    _check_dem(target)
    return target


def _user_dem(path: Path, out_dir: Path, log: Path) -> Path:
    if not path.is_file():
        raise StepError(f"dem_file not found: {path}")
    if path.with_name(path.name + ".xml").is_file():
        _check_dem(path)
        return path
    target = out_dir / (path.stem + ".dem")
    run(["gdal_translate", "-of", "ENVI", "-ot", "Float32", str(path), str(target)], log)
    run(["gdal2isce_xml.py", "-i", str(target)], log)
    _check_dem(target)
    return target


def _check_dem(path: Path) -> None:
    """A failed tile download leaves zeros: refuse such a DEM."""
    from osgeo import gdal
    ds = gdal.Open(str(path))
    if ds is None:
        raise StepError(f"Cannot open DEM {path}")
    a = ds.GetRasterBand(1).ReadAsArray().astype(float)
    good = np.isfinite(a) & (a != 0) & (a > -500)
    share = good.mean()
    if share < 0.98 or np.ptp(a[good]) < 50:
        raise StepError(f"DEM {path} looks wrong: {share:.1%} valid cells, range "
                        f"{np.nanmin(a):.0f}..{np.nanmax(a):.0f} m (a blocked tile download leaves zeros). "
                        "Delete it and rerun, or set isce.dem_file.")
    print(f"  DEM ok: {ds.RasterXSize} x {ds.RasterYSize}, {np.min(a[good]):.0f}..{np.max(a[good]):.0f} m")


def aux_dir(cfg: Config) -> Path:
    """ISCE2 needs an aux-calibration folder; for IPF >= 2.36 (all data since 2016) it may be empty."""
    d = cfg.work_dir / "aux_cal"
    d.mkdir(parents=True, exist_ok=True)
    return d


def run_step(cfg: Config, key: str) -> None:
    orbits(cfg, key)
    dem(cfg)
    aux_dir(cfg)
    _ = tool_env()
