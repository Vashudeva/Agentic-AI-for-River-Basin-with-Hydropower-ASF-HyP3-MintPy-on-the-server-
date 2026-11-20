"""Step 'products': geocoded MintPy grids -> point products in the agent's format.

For every tropospheric variant it writes, in <export_dir>/points/<variant>/:

  os_native_<product>.gpkg       what the open-source chain finds on its own: pixels with
                                 temporal coherence >= export.min_tcoh, thinned to one per
                                 export.thin_cell_m cell (most coherent kept)
  os_at_sarscape_<product>.gpkg  the open-source time series at every SARscape point
                                 (same geometry, attributes and elevations as the SARscape
                                 file; displacement, velocity and coherence replaced)
  sarscape_matched_<product>.gpkg  the SARscape points that received a value, unchanged,
                                 restricted to the common dates - the like-for-like partner
  offset_scan.csv                geolocation check: spread of SARscape - open-source velocity
                                 for shifted sampling positions
  export_summary.json

products: ascending, descending (LOS) and vertical, east (decomposition of the two LOS
grids with d_N = 0, as in SARscape; each LOS series is interpolated to the union of the
dates within their common span). Displacements are in mm, re-zeroed at the first date,
positive towards the satellite (LOS), up (vertical) and east.

It also writes the agent settings files for the comparison runs (<export_dir>/agent_settings).
"""
from __future__ import annotations

import json
import re
import tomllib
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

import h5py
import numpy as np
import pandas as pd

from .config import DECOMPOSED, TRACKS, Config, to_windows

DATE_FIELD = re.compile(r"D_\d{8}$")
METRIC_CRS = "EPSG:32644"          # UTM 44N, as in the agent settings


# --------------------------------------------------------------------------- grids
@dataclass
class Grid:
    key: str
    dates: list[str]
    D: np.ndarray            # (n_dates, rows, cols) mm, re-zeroed at the first date
    tcoh: np.ndarray
    height: np.ndarray
    inc: np.ndarray
    az: np.ndarray
    vel: np.ndarray          # mm/yr
    ps: np.ndarray           # bool
    y0: float
    x0: float
    dy: float
    dx: float
    ref_lalo: tuple[float, float]

    @property
    def shape(self):
        return self.tcoh.shape

    def lat(self, i):
        return self.y0 + (np.asarray(i) + 0.5) * self.dy

    def lon(self, j):
        return self.x0 + (np.asarray(j) + 0.5) * self.dx

    def index(self, lon, lat):
        i = np.floor((np.asarray(lat) - self.y0) / self.dy).astype(int)
        j = np.floor((np.asarray(lon) - self.x0) / self.dx).astype(int)
        return i, j

    def pixel_m(self) -> tuple[float, float]:
        latc = self.y0 + self.dy * self.shape[0] / 2
        return abs(self.dy) * 111_320.0, abs(self.dx) * 111_320.0 * np.cos(np.radians(latc))


def _years(dates: list[str]) -> np.ndarray:
    d = pd.to_datetime(dates, format="%Y%m%d")
    return ((d - d[0]).days.to_numpy(float)) / 365.25


def _first_dataset(f: h5py.File, prefer: str):
    if prefer in f:
        return f[prefer][()]
    for k in f.keys():
        if isinstance(f[k], h5py.Dataset) and f[k].ndim == 2:
            return f[k][()]
    raise KeyError(prefer)


def _projected_epsg(a: dict) -> int | None:
    """EPSG code if the MintPy file is in a projected system (HyP3 products are in UTM)."""
    get = lambda k: (a[k].decode() if isinstance(a.get(k), bytes) else a.get(k))
    epsg = get("EPSG")
    if epsg not in (None, "", "None"):
        e = int(float(epsg))
        return None if e == 4326 else e
    zone = get("UTM_ZONE")
    if zone:
        z = str(zone).strip().upper()
        return (32700 if z.endswith("S") else 32600) + int("".join(ch for ch in z if ch.isdigit()))
    return None


class _ToLatLon:
    """Nearest-neighbour resampling of a projected MintPy grid onto the lat/lon grid of
    config mintpy.geocode_snwe / geocode_step_deg (the grid of the SARscape products)."""

    def __init__(self, cfg: Config, a: dict, epsg: int):
        from pyproj import Transformer
        f = lambda k: float(a[k].decode() if isinstance(a[k], bytes) else a[k])
        S, N, W, E = (float(v) for v in cfg.mintpy["geocode_snwe"])
        step = float(cfg.mintpy.get("geocode_step_deg", 0.000125))
        self.rows, self.cols = int(round((N - S) / step)), int(round((E - W) / step))
        self.y0, self.x0, self.dy, self.dx = N, W, -step, step
        lat = N - (np.arange(self.rows) + 0.5) * step
        lon = W + (np.arange(self.cols) + 0.5) * step
        LON, LAT = np.meshgrid(lon, lat)
        X, Y = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True).transform(LON, LAT)
        sy0, sx0, sdy, sdx = f("Y_FIRST"), f("X_FIRST"), f("Y_STEP"), f("X_STEP")
        length, width = int(f("LENGTH")), int(f("WIDTH"))
        i = np.floor((Y - sy0) / sdy).astype(np.int64)
        j = np.floor((X - sx0) / sdx).astype(np.int64)
        self.ok = (i >= 0) & (i < length) & (j >= 0) & (j < width)
        self.i, self.j = i[self.ok], j[self.ok]
        self.ref = None
        if "REF_Y" in a and "REF_X" in a:
            ry, rx = f("REF_Y"), f("REF_X")
            back = Transformer.from_crs(f"EPSG:{epsg}", "EPSG:4326", always_xy=True)
            rlon, rlat = back.transform(sx0 + (rx + 0.5) * sdx, sy0 + (ry + 0.5) * sdy)
            self.ref = (float(rlat), float(rlon))

    def __call__(self, arr: np.ndarray, fill=np.nan) -> np.ndarray:
        if arr.ndim == 2:
            out = np.full((self.rows, self.cols), fill, dtype=np.float32)
            out[self.ok] = arr[self.i, self.j]
        else:
            out = np.full((arr.shape[0], self.rows, self.cols), fill, dtype=np.float32)
            out[:, self.ok] = arr[:, self.i, self.j]
        return out


def read_grid(cfg: Config, key: str, variant: str) -> Grid:
    """Geocoded MintPy files of one variant. Grids in a projected system (HyP3, UTM) are
    resampled to the lat/lon grid of the SARscape products, so everything downstream is
    the same for both processing routes."""
    d = cfg.geo_dir(key, variant)
    with h5py.File(d / "geo_timeseries.h5", "r") as f:
        D = f["timeseries"][()].astype(np.float32) * 1000.0
        dates = [x.decode() if isinstance(x, bytes) else str(x) for x in f["date"][()]]
        a = dict(f.attrs)
    order = np.argsort(dates)
    D, dates = D[order], [dates[k] for k in order]
    with h5py.File(d / "geo_temporalCoherence.h5", "r") as f:
        tcoh = _first_dataset(f, "temporalCoherence").astype(np.float32)
    with h5py.File(d / "geo_geometryRadar.h5", "r") as f:
        height = f["height"][()].astype(np.float32)
        inc = f["incidenceAngle"][()].astype(np.float32)
        az = f["azimuthAngle"][()].astype(np.float32)
    vel = None
    if (d / "geo_velocity.h5").is_file():
        with h5py.File(d / "geo_velocity.h5", "r") as f:
            if "velocity" in f:
                vel = f["velocity"][()].astype(np.float32) * 1000.0
    ps = None
    if (d / "geo_maskPS.h5").is_file():
        with h5py.File(d / "geo_maskPS.h5", "r") as f:
            ps = np.nan_to_num(_first_dataset(f, "mask")).astype(np.float32)

    attr = lambda k: float(a[k].decode() if isinstance(a[k], bytes) else a[k])
    epsg = _projected_epsg(a)
    if epsg:
        rs = _ToLatLon(cfg, a, epsg)
        D, tcoh, height, inc, az = rs(D), rs(tcoh, 0.0), rs(height), rs(inc), rs(az)
        vel = rs(vel) if vel is not None else None
        ps = rs(ps, 0.0) if ps is not None else None
        y0, x0, dy, dx = rs.y0, rs.x0, rs.dy, rs.dx
        ref = rs.ref or (np.nan, np.nan)
    else:
        y0, x0, dy, dx = attr("Y_FIRST"), attr("X_FIRST"), attr("Y_STEP"), attr("X_STEP")
        ref = (attr("REF_LAT"), attr("REF_LON")) if "REF_LAT" in a else (np.nan, np.nan)
    D -= D[0]
    if vel is None:
        vel = _ols_grid(D, _years(dates))
    ps = (ps > 0) if ps is not None else np.zeros(tcoh.shape, bool)
    tcoh = np.nan_to_num(tcoh)
    tcoh[~np.isfinite(D).all(axis=0)] = 0
    return Grid(key, dates, D, tcoh, height, inc, az, vel, ps, y0, x0, dy, dx, ref)


def _ols_grid(D: np.ndarray, t: np.ndarray) -> np.ndarray:
    tc = t - t.mean()
    return (np.tensordot(tc, D - np.nanmean(D, axis=0), axes=(0, 0)) / (tc @ tc)).astype(np.float32)


def ols(series: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Velocity per row of series (n_points, n_dates)."""
    tc = t - t.mean()
    return (series - series.mean(axis=1, keepdims=True)) @ tc / (tc @ tc)


def common_grid(a: Grid, b: Grid) -> tuple[Grid, Grid]:
    """Crop two grids with the same step to their common extent."""
    if not (np.isclose(a.dy, b.dy) and np.isclose(a.dx, b.dx)):
        raise ValueError("ascending and descending grids have different pixel sizes; geocode both "
                         "with the same mintpy.geocode_step_deg")
    di = int(round((b.y0 - a.y0) / a.dy))
    dj = int(round((b.x0 - a.x0) / a.dx))
    ia0, ja0 = max(0, di), max(0, dj)
    ib0, jb0 = max(0, -di), max(0, -dj)
    n = min(a.shape[0] - ia0, b.shape[0] - ib0)
    m = min(a.shape[1] - ja0, b.shape[1] - jb0)

    def crop(g: Grid, i0, j0):
        s = (slice(i0, i0 + n), slice(j0, j0 + m))
        return Grid(g.key, g.dates, g.D[(slice(None),) + s], g.tcoh[s], g.height[s], g.inc[s], g.az[s],
                    g.vel[s], g.ps[s], g.y0 + i0 * g.dy, g.x0 + j0 * g.dx, g.dy, g.dx, g.ref_lalo)

    return crop(a, ia0, ja0), crop(b, ib0, jb0)


# --------------------------------------------------------------------------- decomposition
def unit_vectors(inc: np.ndarray, az: np.ndarray):
    """(east, up) components of the LOS unit vector (ground -> satellite), MintPy convention."""
    try:
        from mintpy.utils.utils0 import get_unit_vector4component_of_interest as uv
        e, _, u = uv(inc, az, comp="enu2los")
    except ImportError:
        e = -np.sin(np.radians(inc)) * np.sin(np.radians(az))
        u = np.cos(np.radians(inc))
    return np.asarray(e, float), np.asarray(u, float)


def union_dates(a: list[str], b: list[str]) -> list[str]:
    lo, hi = max(a[0], b[0]), min(a[-1], b[-1])
    return sorted({d for d in a + b if lo <= d <= hi})


def interp_matrix(src: list[str], dst: list[str]) -> np.ndarray:
    ts = pd.to_datetime(src, format="%Y%m%d").to_numpy("datetime64[D]").astype(float)
    td = pd.to_datetime(dst, format="%Y%m%d").to_numpy("datetime64[D]").astype(float)
    eye = np.eye(len(src))
    return np.column_stack([np.interp(td, ts, eye[:, c]) for c in range(len(src))])


def decompose(A: Grid, B: Grid, ii: np.ndarray, jj: np.ndarray):
    """Vertical and east histories (n_pix, n_union) at pixels (ii, jj) of the common grid."""
    du = union_dates(A.dates, B.dates)
    WA, WB = interp_matrix(A.dates, du), interp_matrix(B.dates, du)
    dA = WA @ A.D[:, ii, jj]                  # (n_union, n_pix)
    dB = WB @ B.D[:, ii, jj]
    eA, uA = unit_vectors(A.inc[ii, jj], A.az[ii, jj])
    eB, uB = unit_vectors(B.inc[ii, jj], B.az[ii, jj])
    det = eA * uB - eB * uA
    E = (dA * uB - dB * uA) / det
    U = (eA * dB - eB * dA) / det
    E, U = E - E[0], U - U[0]
    return du, U.T, E.T


def decompose_velocity(A: Grid, B: Grid):
    eA, uA = unit_vectors(A.inc, A.az)
    eB, uB = unit_vectors(B.inc, B.az)
    det = eA * uB - eB * uA
    with np.errstate(invalid="ignore", divide="ignore"):
        return (eA * B.vel - eB * A.vel) / det, (A.vel * uB - B.vel * uA) / det     # up, east


# --------------------------------------------------------------------------- geometry helpers
def to_metric(lon, lat):
    from pyproj import Transformer
    tr = Transformer.from_crs("EPSG:4326", METRIC_CRS, always_xy=True)
    return tr.transform(np.asarray(lon, float), np.asarray(lat, float))


def thin(lon, lat, score, cell_m: float, max_points: int) -> np.ndarray:
    """Indices of the best-scoring point per cell, at most max_points (highest scores)."""
    x, y = to_metric(lon, lat)
    cell = pd.DataFrame({"cx": np.floor(np.asarray(x) / cell_m).astype(np.int64),
                         "cy": np.floor(np.asarray(y) / cell_m).astype(np.int64),
                         "s": score, "k": np.arange(len(score))})
    best = cell.sort_values("s", ascending=False).drop_duplicates(["cx", "cy"])
    return np.sort(best["k"].to_numpy()[:max_points])


def search_offsets(radius_px_y: float, radius_px_x: float):
    ry, rx = int(np.ceil(radius_px_y)), int(np.ceil(radius_px_x))
    off = [(di, dj, np.hypot(di / max(radius_px_y, 1e-9), dj / max(radius_px_x, 1e-9)))
           for di in range(-ry, ry + 1) for dj in range(-rx, rx + 1)]
    return [(di, dj) for di, dj, r in sorted(off, key=lambda o: o[2]) if r <= 1.0 + 1e-9]


def pick_pixels(g: Grid, quality: np.ndarray, lon, lat, min_q: float, radius_m: float):
    """Pixel for each position: the containing pixel if quality >= min_q, else the best
    pixel within radius_m. Returns i, j, ok, distance (m)."""
    i0, j0 = g.index(lon, lat)
    rows, cols = g.shape
    n = len(i0)

    def q_at(i, j):
        inside = (i >= 0) & (i < rows) & (j >= 0) & (j < cols)
        q = np.full(n, -np.inf)
        q[inside] = quality[i[inside], j[inside]]
        return q

    q0 = q_at(i0, j0)
    fixed = q0 >= min_q                         # the containing pixel qualifies: keep it
    bi, bj = np.where(fixed, i0, -1), np.where(fixed, j0, -1)
    bq = np.where(fixed, q0, -np.inf)
    py, px = g.pixel_m()
    for di, dj in search_offsets(radius_m / py, radius_m / px):
        if di == 0 and dj == 0:
            continue
        i, j = i0 + di, j0 + dj
        q = q_at(i, j)
        take = ~fixed & (q >= min_q) & (q > bq)
        bi[take], bj[take], bq[take] = i[take], j[take], q[take]
    ok = bi >= 0
    dist = np.full(n, np.nan)
    if ok.any():
        x0, y0 = to_metric(np.asarray(lon)[ok], np.asarray(lat)[ok])
        x1, y1 = to_metric(g.lon(bj[ok]), g.lat(bi[ok]))
        dist[ok] = np.hypot(np.asarray(x1) - x0, np.asarray(y1) - y0)
    return bi, bj, ok, dist


# --------------------------------------------------------------------------- writing
def write_points(path: Path, frame: pd.DataFrame, lon, lat, geometry=None) -> None:
    import geopandas as gpd
    geom = geometry if geometry is not None else gpd.points_from_xy(lon, lat)
    g = gpd.GeoDataFrame(frame.reset_index(drop=True), geometry=geom, crs="EPSG:4326")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    g.to_file(path, driver="GPKG", layer=path.stem, engine="pyogrio")


def series_frame(dates: list[str], S: np.ndarray) -> pd.DataFrame:
    return pd.DataFrame(np.round(S, 3).astype(np.float64), columns=[f"D_{d}" for d in dates])


def native_los(cfg: Config, g: Grid) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    e = cfg.export
    good = (g.tcoh >= float(e.get("min_tcoh", 0.7))) & np.isfinite(g.height) & np.isfinite(g.vel)
    ii, jj = np.nonzero(good)
    lon, lat = g.lon(jj), g.lat(ii)
    keep = thin(lon, lat, g.tcoh[ii, jj], float(e.get("thin_cell_m", 30)), int(e.get("max_points", 60000)))
    ii, jj, lon, lat = ii[keep], jj[keep], lon[keep], lat[keep]
    f = series_frame(g.dates, g.D[:, ii, jj].T)
    f["velocity"] = np.round(g.vel[ii, jj], 3)
    f["coherence"] = np.round(g.tcoh[ii, jj], 3)
    f["Z"] = np.round(g.height[ii, jj], 2)
    f["scatterer"] = np.where(g.ps[ii, jj], "PS", "DS")
    f["ILOS"] = np.round(g.inc[ii, jj], 3)
    f["AZ_mintpy"] = np.round(g.az[ii, jj], 3)
    f["lon"], f["lat"] = lon, lat
    return f, lon, lat


def native_decomposed(cfg: Config, A: Grid, B: Grid):
    e = cfg.export
    q = np.minimum(A.tcoh, B.tcoh)
    good = (q >= float(e.get("min_tcoh", 0.7))) & np.isfinite(A.height)
    ii, jj = np.nonzero(good)
    lon, lat = A.lon(jj), A.lat(ii)
    keep = thin(lon, lat, q[ii, jj], float(e.get("thin_cell_m", 30)), int(e.get("max_points", 60000)))
    ii, jj, lon, lat = ii[keep], jj[keep], lon[keep], lat[keep]
    du, U, E = decompose(A, B, ii, jj)
    t = _years(du)
    out = {}
    for name, S in (("vertical", U), ("east", E)):
        f = series_frame(du, S)
        f["velocity"] = np.round(ols(S, t), 3)
        f["coherence"] = np.round(q[ii, jj], 3)
        f["h"] = np.round(A.height[ii, jj], 2)
        f["scatterer"] = [f"{'PS' if a else 'DS'}-{'PS' if b else 'DS'}" for a, b in zip(A.ps[ii, jj], B.ps[ii, jj])]
        f["ascILOS"], f["dscILOS"] = np.round(A.inc[ii, jj], 3), np.round(B.inc[ii, jj], 3)
        f["lon"], f["lat"] = lon, lat
        out[name] = f
    return out, lon, lat


def read_sarscape(path: Path):
    import pyogrio
    fr = pyogrio.read_dataframe(path)
    dcols = sorted(c for c in fr.columns if DATE_FIELD.fullmatch(c))
    if "lon" in fr.columns and "lat" in fr.columns:
        lon = pd.to_numeric(fr["lon"], errors="coerce").to_numpy(float)
        lat = pd.to_numeric(fr["lat"], errors="coerce").to_numpy(float)
        src = "lon/lat fields"
    else:
        g = fr.geometry.to_crs(4326) if fr.crs and fr.crs.to_epsg() != 4326 else fr.geometry
        lon, lat, src = g.x.to_numpy(float), g.y.to_numpy(float), "geometry"
    return fr, dcols, lon, lat, src


def at_sarscape(cfg: Config, key: str, fr, dcols, lon, lat, dates_os, S_os, ok, extra: dict):
    """Two frames on identical points/dates: SARscape (matched) and open-source values."""
    sar_dates = [c[2:] for c in dcols]
    common = [d for d in sar_dates if d in set(dates_os)]
    if len(common) < 3:
        raise ValueError(f"{key}: fewer than three common dates between SARscape and the open-source run")
    t = _years(common)
    k_os = [dates_os.index(d) for d in common]
    rows = np.flatnonzero(ok)
    base = fr.iloc[rows].copy()
    other = [c for c in fr.columns if c not in dcols and c != "geometry"]
    sar = base[other].reset_index(drop=True)
    Dsar = base[[f"D_{d}" for d in common]].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    Dsar = Dsar - Dsar[:, [0]]
    matched = pd.concat([sar, series_frame(common, Dsar)], axis=1)
    matched["sar_row"] = rows

    Dos = S_os[rows][:, k_os]
    Dos = Dos - Dos[:, [0]]
    os_frame = pd.concat([sar.copy(), series_frame(common, Dos)], axis=1)
    os_frame["velocity"] = np.round(ols(Dos, t), 3)
    for name, values in extra.items():
        os_frame[name] = values[rows] if isinstance(values, np.ndarray) else values
    if "coherence" in extra:
        os_frame["coherence"] = extra["coherence"][rows]
    os_frame["sar_row"] = rows
    geometry = base.geometry.reset_index(drop=True)
    return matched, os_frame, geometry, common


# --------------------------------------------------------------------------- offset scan
def offset_scan(cfg: Config, name: str, vel_grid: np.ndarray, q_grid: np.ndarray, g: Grid,
                lon, lat, v_sar, positions: str) -> pd.DataFrame:
    e = cfg.export
    half, step = float(e.get("offset_scan_m", 45)), float(e.get("offset_step_m", 5))
    thr = float(e.get("min_tcoh", 0.7))
    py, px = g.pixel_m()
    rows = []
    shifts = np.arange(-half, half + 1e-6, step)
    for sy in shifts:
        for sx in shifts:
            i, j = g.index(np.asarray(lon) + sx / px * abs(g.dx), np.asarray(lat) + sy / py * abs(g.dy))
            inside = (i >= 0) & (i < g.shape[0]) & (j >= 0) & (j < g.shape[1])
            v = np.full(len(i), np.nan)
            q = np.zeros(len(i))
            v[inside] = vel_grid[i[inside], j[inside]]
            q[inside] = q_grid[i[inside], j[inside]]
            d = v - v_sar
            use = (q >= thr) & np.isfinite(d)
            if use.sum() < 50:
                continue
            med = np.median(d[use])
            rows.append({"product": name, "positions": positions, "east_m": sx, "north_m": sy, "points": int(use.sum()),
                         "median_difference": med, "nmad": 1.4826 * np.median(np.abs(d[use] - med))})
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- agent settings
def _toml_value(v) -> str:
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, (int, float)):
        return repr(v)
    if isinstance(v, (list, tuple)):
        return "[" + ", ".join(_toml_value(x) for x in v) + "]"
    s = str(v)
    return f"'{s}'" if "'" not in s else json.dumps(s)


def write_agent_settings(cfg: Config, path: Path, products: dict, work_root: Path, note: str) -> None:
    src = cfg.agent_dir / "settings.toml" if cfg.agent_dir else None
    raw = tomllib.loads(src.read_text()) if src and src.is_file() else {}
    sites = raw.get("paths", {}).get("sites", "sites/sites.geojson")
    masks = raw.get("paths", {}).get("masks", "sites/masks.geojson")
    absolute = lambda p: to_windows(cfg.agent_dir / str(p).replace("\\", "/")) if cfg.agent_dir and not Path(
        str(p)).is_absolute() else str(p)
    lines = [f"# {note}", "# Written by kgproc (open-source comparison). Use with:  set KALI_SETTINGS=<this file>", "",
             "[paths]", f"data_root = {_toml_value(to_windows(cfg.export_dir))}",
             f"work_root = {_toml_value(to_windows(work_root))}",
             f"sites = {_toml_value(absolute(sites))}", f"masks = {_toml_value(absolute(masks))}", ""]
    for key, p in products.items():
        lines.append(f"[products.{key}]")
        for k, v in p.items():
            lines.append(f"{k} = {_toml_value(v)}")
        lines.append("")
    for section in ("site", "screening", "elevation_check", "zones"):
        if section in raw:
            lines.append(f"[{section}]")
            lines += [f"{k} = {_toml_value(v)}" for k, v in raw[section].items()]
            lines.append("")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines))


def agent_product_defaults(cfg: Config) -> dict:
    """mask, shift_deg and elevation_field per product, from the agent's settings.toml."""
    src = cfg.agent_dir / "settings.toml" if cfg.agent_dir else None
    raw = tomllib.loads(src.read_text()) if src and src.is_file() else {}
    out = {}
    for key in TRACKS + DECOMPOSED:
        p = raw.get("products", {}).get(key, {})
        out[key] = {"kind": p.get("kind", "los" if key in TRACKS else key),
                    "label": p.get("label", key),
                    "elevation_field": p.get("elevation_field", "Z" if key in TRACKS else "h"),
                    "mask": p.get("mask", ""), "shift_deg": p.get("shift_deg", [0.0, 0.0]),
                    "file": p.get("file", "")}
    return out, raw


# --------------------------------------------------------------------------- main
def run_step(cfg: Config, variant: str) -> dict:
    out = cfg.points_dir(variant)
    out.mkdir(parents=True, exist_ok=True)
    e = cfg.export
    smin, srad = float(e.get("sample_min_tcoh", 0.5)), float(e.get("sample_search_m", 15))
    A = read_grid(cfg, "ascending", variant)
    B = read_grid(cfg, "descending", variant)
    if not np.allclose(A.ref_lalo, B.ref_lalo, atol=max(2 * abs(A.dy), 0.0005)):
        print(f"  WARNING: reference points differ (asc {A.ref_lalo}, dsc {B.ref_lalo}); "
              "the decomposition needs one common reference")
    summary = {"variant": variant, "created": datetime.now().isoformat(timespec="seconds"),
               "reference_lalo": {"ascending": A.ref_lalo, "descending": B.ref_lalo}, "products": {}}
    scans = []
    defaults, _ = agent_product_defaults(cfg)

    # ---- line of sight
    for g in (A, B):
        key = g.key
        f, lon, lat = native_los(cfg, g)
        write_points(out / f"os_native_{key}.gpkg", f, lon, lat)
        fr, dcols, slon, slat, src = read_sarscape(cfg.sarscape_products[key])
        i, j, ok, dist = pick_pixels(g, g.tcoh, slon, slat, smin, srad)
        S_os = np.full((len(fr), len(g.dates)), np.nan)
        S_os[ok] = g.D[:, i[ok], j[ok]].T
        extra = {"coherence": np.where(ok, g.tcoh[np.clip(i, 0, None), np.clip(j, 0, None)], np.nan),
                 "os_dist_m": dist, "os_height": np.where(ok, g.height[np.clip(i, 0, None), np.clip(j, 0, None)], np.nan),
                 "os_ps": np.where(ok, g.ps[np.clip(i, 0, None), np.clip(j, 0, None)], False).astype(int)}
        matched, osf, geom, common = at_sarscape(cfg, key, fr, dcols, slon, slat, g.dates, S_os, ok, extra)
        write_points(out / f"sarscape_matched_{key}.gpkg", matched, None, None, geometry=geom)
        write_points(out / f"os_at_sarscape_{key}.gpkg", osf, None, None, geometry=geom)
        v_sar = pd.to_numeric(fr["velocity"], errors="coerce").to_numpy(float)
        scans.append(offset_scan(cfg, key, g.vel, g.tcoh, g, slon, slat, v_sar, src))
        summary["products"][key] = {"native_points": len(f), "sarscape_points": len(fr), "matched": int(ok.sum()),
                                    "common_dates": len(common), "open_source_dates": len(g.dates),
                                    "sarscape_dates": len(dcols), "sampling_positions": src,
                                    "median_sampling_distance_m": float(np.nanmedian(dist)) if ok.any() else None}
        print(f"  {key}: native {len(f)} points; {int(ok.sum())}/{len(fr)} SARscape points matched; "
              f"{len(common)} common dates")

    # ---- decomposition on the common grid
    A2, B2 = common_grid(A, B)
    out_dec, lon, lat = native_decomposed(cfg, A2, B2)
    for name in DECOMPOSED:
        write_points(out / f"os_native_{name}.gpkg", out_dec[name], lon, lat)
    vU, vE = decompose_velocity(A2, B2)
    q2 = np.minimum(A2.tcoh, B2.tcoh)
    for name in DECOMPOSED:
        if name not in cfg.sarscape_products:
            continue
        fr, dcols, slon, slat, src = read_sarscape(cfg.sarscape_products[name])
        i, j, ok, dist = pick_pixels(A2, q2, slon, slat, smin, srad)
        du, U, E = decompose(A2, B2, i[ok], j[ok])
        S_os = np.full((len(fr), len(du)), np.nan)
        S_os[ok] = U if name == "vertical" else E
        cl = lambda arr: np.where(ok, arr[np.clip(i, 0, None), np.clip(j, 0, None)], np.nan)
        extra = {"coherence": cl(q2), "os_dist_m": dist, "os_height": cl(A2.height)}
        matched, osf, geom, common = at_sarscape(cfg, name, fr, dcols, slon, slat, du, S_os, ok, extra)
        write_points(out / f"sarscape_matched_{name}.gpkg", matched, None, None, geometry=geom)
        write_points(out / f"os_at_sarscape_{name}.gpkg", osf, None, None, geometry=geom)
        v_sar = pd.to_numeric(fr["velocity"], errors="coerce").to_numpy(float)
        vgrid = vU if name == "vertical" else vE
        scans.append(offset_scan(cfg, name, vgrid, q2, A2, slon, slat, v_sar, src))
        if src != "geometry":                     # also scan the shapefile geometry positions
            gx, gy = fr.geometry.x.to_numpy(float), fr.geometry.y.to_numpy(float)
            if np.nanmax(np.hypot(gx - slon, gy - slat)) > 1e-6:
                scans.append(offset_scan(cfg, name, vgrid, q2, A2, gx, gy, v_sar, "geometry"))
        summary["products"][name] = {"native_points": len(out_dec[name]), "sarscape_points": len(fr),
                                     "matched": int(ok.sum()), "common_dates": len(common),
                                     "open_source_dates": len(du), "sarscape_dates": len(dcols),
                                     "sampling_positions": src,
                                     "median_sampling_distance_m": float(np.nanmedian(dist)) if ok.any() else None}
        print(f"  {name}: native {len(out_dec[name])} points; {int(ok.sum())}/{len(fr)} SARscape points matched; "
              f"{len(common)} common dates")

    scan = pd.concat(scans, ignore_index=True) if scans else pd.DataFrame()
    scan.to_csv(out / "offset_scan.csv", index=False)
    if len(scan):
        best = scan.loc[scan.groupby(["product", "positions"])["nmad"].idxmin()]
        summary["offset_scan_best"] = best.to_dict(orient="records")
    (out / "export_summary.json").write_text(json.dumps(summary, indent=1, default=float))
    settings_files(cfg, variant, summary)
    return summary


def settings_files(cfg: Config, variant: str, summary: dict) -> None:
    """Agent settings for the comparison runs (one set of four products each)."""
    defaults, raw = agent_product_defaults(cfg)
    sdir = cfg.export_dir / "agent_settings"
    runs = cfg.export_dir / "agent_runs"
    pts = cfg.points_dir(variant)
    data_root = str(raw.get("paths", {}).get("data_root") or to_windows(cfg.sarscape_root))

    def sarscape_file(key) -> str:
        f = defaults[key]["file"]
        if not f:
            return to_windows(cfg.sarscape_products[key])
        if re.match(r"^[A-Za-z]:[\\/]", f) or f.startswith("/"):
            return f
        sep = "\\" if re.match(r"^[A-Za-z]:", data_root) else "/"
        return data_root.rstrip("\\/") + sep + f.replace("/", sep).replace("\\", sep)

    def entry(key, file, shift, n, label_suffix, elevation=None):
        d = defaults[key]
        base = d["label"] if "SARscape" in label_suffix else d["label"].replace(", E-SBAS", "")
        e = {"file": file if isinstance(file, str) else to_windows(file), "kind": d["kind"],
             "label": f"{base} {label_suffix}", "elevation_field": elevation or d["elevation_field"]}
        if d["mask"]:
            e["mask"] = d["mask"]
        if n is not None:
            e["expected_points"] = int(n)
        if shift and any(shift):
            e["shift_deg"] = list(shift)
        return e

    keys = [k for k in TRACKS + DECOMPOSED if k in summary["products"]]
    # 1. SARscape as delivered (reproduces the manuscript numbers)
    prods = {k: entry(k, sarscape_file(k), defaults[k]["shift_deg"], None, "(SARscape)") for k in keys}
    write_agent_settings(cfg, sdir / "sarscape.toml", prods, runs / "sarscape",
                         "SARscape E-SBAS products as used for the manuscript")
    # 2-4. matched SARscape subset, open-source at the SARscape points, open-source native
    sets = {
        f"sarscape_matched_{variant}": ("sarscape_matched_{}.gpkg", True, "(SARscape, matched points)"),
        f"opensource_at_sarscape_{variant}": ("os_at_sarscape_{}.gpkg", True, f"(open source, {variant}, at SARscape points)"),
        f"opensource_native_{variant}": ("os_native_{}.gpkg", False, f"(open source, {variant}, native points)"),
    }
    for name, (pattern, same_geometry, suffix) in sets.items():
        prods = {}
        for k in keys:
            n = summary["products"][k]["matched" if same_geometry else "native_points"]
            shift = defaults[k]["shift_deg"] if same_geometry else None
            elev = None if same_geometry else ("Z" if k in TRACKS else "h")
            prods[k] = entry(k, pts / pattern.format(k), shift, n, suffix, elev)
        write_agent_settings(cfg, sdir / f"{name}.toml", prods, runs / name, f"{name}: {suffix}")
    print(f"  agent settings in {sdir}")
