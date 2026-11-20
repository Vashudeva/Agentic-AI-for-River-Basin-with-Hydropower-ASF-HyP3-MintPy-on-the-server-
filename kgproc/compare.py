"""Step 'compare': SARscape E-SBAS versus the open-source chain, point by point.

Inputs are the files of the 'products' step (one variant). For each product (ascending,
descending, vertical, east) the SARscape points and the open-source values at the same
points are compared on the same dates, after both are re-referenced in the same way
(compare.reference: "headworks" = median history of the points in zones A-F, as in the
manuscript's double differences; "window" = median of all points).

Outputs in <export_dir>/comparison/<variant>/:
  comparison_report.md          tables and figure links
  velocity_agreement.csv        bias, spread, correlation, slope, shares within 1 / 2 mm/yr
  timeseries_agreement.csv      per-point RMS difference of the histories
  common_mode_<product>.csv     median difference per date (what differs everywhere at once)
  elevation_gradient.csv        Theil-Sen velocity-height slope, far field, both chains
  points_<product>.csv          per-point velocities and differences (for QGIS)
  fig_<product>.png, fig_elevation_gradient.png

Only numpy, pandas, scipy, matplotlib and pyogrio are needed, so this step also runs in
the agent's Windows environment.
"""
from __future__ import annotations

import json
import re
import tomllib
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from .config import DECOMPOSED, TRACKS, Config, local_path

DATE_FIELD = re.compile(r"D_\d{8}$")


# --------------------------------------------------------------------------- reading
def read_points(path: Path):
    import pyogrio
    fr = pyogrio.read_dataframe(path)
    dcols = sorted(c for c in fr.columns if DATE_FIELD.fullmatch(c))
    D = fr[dcols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    D = D - D[:, [0]]
    return fr, dcols, D


def years(dcols: list[str]) -> np.ndarray:
    d = pd.to_datetime([c[2:] for c in dcols], format="%Y%m%d")
    return (d - d[0]).days.to_numpy(float) / 365.25


def ols(D: np.ndarray, t: np.ndarray) -> np.ndarray:
    tc = t - t.mean()
    return (D - np.nanmean(D, axis=1, keepdims=True)) @ tc / (tc @ tc)


def agent_settings(cfg: Config) -> dict:
    f = cfg.agent_dir / "settings.toml" if cfg.agent_dir else None
    return tomllib.loads(f.read_text()) if f and f.is_file() else {}


def zones(cfg: Config, raw: dict):
    """Headworks zones from the agent's site register (shapely polygons)."""
    import shapely
    rel = str(raw.get("paths", {}).get("sites", "sites/sites.geojson"))
    p = local_path(rel)
    if not p.is_absolute():
        p = cfg.agent_dir / rel.replace("\\", "/")
    if not p.is_file():
        return []
    fc = json.loads(p.read_text(encoding="utf-8"))
    group = raw.get("zones", {}).get("group", "headworks")
    out = []
    for f in fc.get("features", []):
        pr = f.get("properties") or {}
        if pr.get("group") == group and f.get("geometry"):
            out.append((pr.get("name"), shapely.geometry.shape(f["geometry"])))
    return out


def positions(fr, shift) -> tuple[np.ndarray, np.ndarray]:
    """Point positions as the agent sees them (geometry + shift_deg)."""
    x, y = fr.geometry.x.to_numpy(float), fr.geometry.y.to_numpy(float)
    return x + float(shift[0]), y + float(shift[1])


def in_zones(lon, lat, zs) -> np.ndarray:
    import shapely
    if not zs:
        return np.zeros(len(lon), bool)
    pts = shapely.points(lon, lat)
    inside = np.zeros(len(lon), bool)
    for _, poly in zs:
        inside |= shapely.contains(poly, pts)
    return inside


def distance_m(lon, lat, lon0, lat0) -> np.ndarray:
    k = 111_320.0
    return np.hypot((np.asarray(lon) - lon0) * k * np.cos(np.radians(lat0)), (np.asarray(lat) - lat0) * k)


# --------------------------------------------------------------------------- statistics
def agreement(v_os: np.ndarray, v_sar: np.ndarray) -> dict:
    ok = np.isfinite(v_os) & np.isfinite(v_sar)
    a, b = v_os[ok], v_sar[ok]
    d = a - b
    med = float(np.median(d))
    cov = np.cov(np.vstack([b, a]))
    evals, evecs = np.linalg.eigh(cov)
    major = evecs[:, np.argmax(evals)]
    return {"points": int(ok.sum()),
            "median_difference_mm_yr": med,
            "nmad_mm_yr": float(1.4826 * np.median(np.abs(d - med))),
            "rmse_mm_yr": float(np.sqrt(np.mean(d ** 2))),
            "pearson_r": float(stats.pearsonr(a, b)[0]) if ok.sum() > 2 else np.nan,
            "spearman_rho": float(stats.spearmanr(a, b)[0]) if ok.sum() > 2 else np.nan,
            "orthogonal_slope": float(major[1] / major[0]) if major[0] != 0 else np.nan,
            "share_within_1_mm_yr": float(np.mean(np.abs(d - med) <= 1.0)),
            "share_within_2_mm_yr": float(np.mean(np.abs(d - med) <= 2.0))}


def gradient(v: np.ndarray, h: np.ndarray, far: np.ndarray, n: int, seed: int) -> dict:
    idx = np.flatnonzero(far & np.isfinite(v) & np.isfinite(h))
    if len(idx) < 20:
        return {"points": int(len(idx))}
    rng = np.random.default_rng(seed)
    idx = rng.choice(idx, min(n, len(idx)), replace=False)
    hk = (h[idx] - np.median(h[idx])) / 1000.0
    s, i, lo, hi = stats.theilslopes(v[idx], hk, 0.95)
    rho = stats.spearmanr(h[idx], v[idx])[0]
    return {"points": int(len(idx)), "slope_mm_yr_per_km": float(s), "ci95_low": float(lo),
            "ci95_high": float(hi), "spearman_rho": float(rho)}


# --------------------------------------------------------------------------- figures
def figure_product(name: str, v_sar, v_os, lon, lat, common, path: Path, hw) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    ok = np.isfinite(v_sar) & np.isfinite(v_os)
    d = v_os - v_sar
    med = np.nanmedian(d[ok])
    nmad = 1.4826 * np.nanmedian(np.abs(d[ok] - med))
    fig, ax = plt.subplots(2, 2, figsize=(11, 9.5))
    lim = np.nanpercentile(np.abs(np.r_[v_sar[ok], v_os[ok]]), 99.5)
    a = ax[0, 0]
    a.hexbin(v_sar[ok], v_os[ok], gridsize=70, bins="log", cmap="Greys", extent=(-lim, lim, -lim, lim), mincnt=1)
    a.plot([-lim, lim], [-lim, lim], color="tab:red", lw=1)
    a.set(xlabel="SARscape E-SBAS velocity (mm/yr)", ylabel="Open-source velocity (mm/yr)",
          title=f"{name}: {ok.sum():,} points", xlim=(-lim, lim), ylim=(-lim, lim))
    a.set_aspect("equal")
    a = ax[0, 1]
    mean = (v_sar + v_os) / 2
    a.scatter(mean[ok], d[ok], s=2, color="0.3", alpha=0.3, lw=0)
    for y, ls in ((med, "-"), (med + 1.96 * nmad, "--"), (med - 1.96 * nmad, "--")):
        a.axhline(y, color="tab:red", ls=ls, lw=1)
    a.set(xlabel="Mean of the two (mm/yr)", ylabel="Open source - SARscape (mm/yr)",
          title=f"Bland-Altman: median {med:+.2f}, NMAD {nmad:.2f} mm/yr")
    a = ax[1, 0]
    sc = a.scatter(lon[ok], lat[ok], c=np.clip(d[ok] - med, -5, 5), s=3, cmap="RdBu_r", vmin=-5, vmax=5, lw=0)
    a.plot(*hw, marker="o", mfc="none", mec="k", ms=10)
    a.set(xlabel="Longitude", ylabel="Latitude", title="Difference minus its median (mm/yr)")
    a.set_aspect(1 / np.cos(np.radians(np.nanmean(lat))))
    fig.colorbar(sc, ax=a, shrink=0.8)
    a = ax[1, 1]
    dates = pd.to_datetime(common["date"])
    a.plot(dates, common["median_difference_mm"], color="k", lw=1, label="median difference")
    a.fill_between(dates, common["median_difference_mm"] - common["nmad_mm"],
                   common["median_difference_mm"] + common["nmad_mm"], color="0.8", label="± NMAD")
    a.axhline(0, color="0.5", lw=0.8)
    a.set(ylabel="Open source - SARscape (mm)", title="Difference of the histories, per date")
    a.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def figure_gradient(rows: list[tuple], path: Path, n: int = 3000, seed: int = 0) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, len(rows), figsize=(4.2 * len(rows), 4), squeeze=False)
    for a, (name, h, vs, vo, far) in zip(axes[0], rows):
        idx = np.flatnonzero(far & np.isfinite(vs) & np.isfinite(vo) & np.isfinite(h))
        idx = np.random.default_rng(seed).choice(idx, min(n, len(idx)), replace=False) if len(idx) else idx
        ok = np.zeros(len(h), bool)
        ok[idx] = True                       # Theil-Sen is O(n^2): same sample size as the table
        a.scatter(h[ok], vs[ok], s=2, color="tab:blue", alpha=0.25, lw=0, label="SARscape")
        a.scatter(h[ok], vo[ok], s=2, color="tab:orange", alpha=0.25, lw=0, label="open source")
        for v, c in ((vs, "tab:blue"), (vo, "tab:orange")):
            if ok.sum() > 20:
                s, i, *_ = stats.theilslopes(v[ok], h[ok] / 1000.0)
                xx = np.array([np.nanmin(h[ok]), np.nanmax(h[ok])])
                a.plot(xx, i + s * xx / 1000.0, color=c, lw=2, label=f"{s:+.1f} mm/yr per km")
        a.set(title=name, xlabel="Elevation (m)", ylabel="Velocity (mm/yr)")
        a.set_ylim(np.nanpercentile(np.r_[vs[ok], vo[ok]], [0.5, 99.5]) if ok.any() else None)
        a.legend(fontsize=7, markerscale=4)
    fig.suptitle("Far field (> far_field_m from the headworks), same points in both chains")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


# --------------------------------------------------------------------------- main
def run_step(cfg: Config, variant: str) -> Path:
    pts = cfg.points_dir(variant)
    out = cfg.report_dir(variant)
    out.mkdir(parents=True, exist_ok=True)
    c = cfg.compare
    raw = agent_settings(cfg)
    zs = zones(cfg, raw)
    hw = (float(c.get("headworks_lon", 83.58004)), float(c.get("headworks_lat", 27.97792)))
    mode = c.get("reference", "headworks")
    far_m, n_far, seed = float(c.get("far_field_m", 1000)), int(c.get("sample_points", 3000)), int(c.get("seed", 0))
    summary = json.loads((pts / "export_summary.json").read_text()) if (pts / "export_summary.json").is_file() else {}

    rows_v, rows_t, rows_g, grad_fig, cover = [], [], [], [], []
    for key in TRACKS + DECOMPOSED:
        fs, fo = pts / f"sarscape_matched_{key}.gpkg", pts / f"os_at_sarscape_{key}.gpkg"
        if not (fs.is_file() and fo.is_file()):
            continue
        S, dcols, Ds = read_points(fs)
        O, dcols_o, Do = read_points(fo)
        if dcols != dcols_o or not np.array_equal(S["sar_row"].to_numpy(), O["sar_row"].to_numpy()):
            raise ValueError(f"{key}: the matched files do not hold the same points and dates")
        p = raw.get("products", {}).get(key, {})
        shift = p.get("shift_deg", [0.0, 0.0])
        elev = p.get("elevation_field", "Z" if key in TRACKS else "h")
        lon, lat = positions(S, shift)
        h = pd.to_numeric(S[elev], errors="coerce").to_numpy(float) if elev in S else np.full(len(S), np.nan)
        t = years(dcols)

        ref_mask = in_zones(lon, lat, zs) if mode == "headworks" else np.ones(len(S), bool)
        if mode == "headworks" and ref_mask.sum() < 5:
            print(f"  {key}: fewer than 5 points in the headworks zones; using the window median")
            ref_mask = np.ones(len(S), bool)
        Rs, Ro = Ds - np.nanmedian(Ds[ref_mask], axis=0), Do - np.nanmedian(Do[ref_mask], axis=0)
        vs, vo = ols(Rs, t), ols(Ro, t)

        a = agreement(vo, vs)
        a.update({"product": key, "variant": variant, "reference": mode, "reference_points": int(ref_mask.sum()),
                  "dates": len(dcols)})
        rows_v.append(a)

        diff = Ro - Rs
        rms = np.sqrt(np.nanmean(diff ** 2, axis=1))
        rows_t.append({"product": key, "points": len(rms), "median_rms_mm": float(np.nanmedian(rms)),
                       "p90_rms_mm": float(np.nanpercentile(rms, 90)),
                       "median_corr": float(np.nanmedian([np.corrcoef(x, y)[0, 1] if np.std(x) > 0 and np.std(y) > 0
                                                          else np.nan for x, y in zip(Rs, Ro)]))})
        med_d = np.nanmedian(diff, axis=0)
        common = pd.DataFrame({"date": pd.to_datetime([d[2:] for d in dcols], format="%Y%m%d").date,
                               "median_difference_mm": med_d,
                               "nmad_mm": 1.4826 * np.nanmedian(np.abs(diff - med_d), axis=0)})
        common.to_csv(out / f"common_mode_{key}.csv", index=False)

        far = distance_m(lon, lat, *hw) > far_m
        gs, go = gradient(vs, h, far, n_far, seed), gradient(vo, h, far, n_far, seed)
        rows_g += [{"product": key, "chain": "SARscape (matched points)", **gs},
                   {"product": key, "chain": "open source (same points)", **go}]
        fn = pts / f"os_native_{key}.gpkg"
        if fn.is_file():
            N, dn, Dn = read_points(fn)
            en = "Z" if key in TRACKS else "h"
            hn = pd.to_numeric(N[en], errors="coerce").to_numpy(float)
            ln, la = N.geometry.x.to_numpy(float), N.geometry.y.to_numpy(float)
            vn = ols(Dn, years(dn))
            rows_g.append({"product": key, "chain": "open source (native points)",
                           **gradient(vn, hn, distance_m(ln, la, *hw) > far_m, n_far, seed)})
        grad_fig.append((key, h, vs, vo, far))

        pd.DataFrame({"sar_row": S["sar_row"], "lon": lon, "lat": lat, "elevation": h,
                      "v_sarscape": vs, "v_open_source": vo, "difference": vo - vs,
                      "rms_history_mm": rms, "os_coherence": O.get("coherence"),
                      "os_dist_m": O.get("os_dist_m")}).to_csv(out / f"points_{key}.csv", index=False)
        figure_product(key, vs, vo, lon, lat, common, out / f"fig_{key}.png", hw)
        s = summary.get("products", {}).get(key, {})
        cover.append({"product": key, "SARscape points": s.get("sarscape_points"), "matched": s.get("matched"),
                      "share matched": (s["matched"] / s["sarscape_points"]) if s.get("sarscape_points") else None,
                      "open-source native points": s.get("native_points"), "common dates": s.get("common_dates"),
                      "SARscape dates": s.get("sarscape_dates"), "open-source dates": s.get("open_source_dates"),
                      "median sampling distance (m)": s.get("median_sampling_distance_m")})
        print(f"  {key}: {a['points']:,} points, median diff {a['median_difference_mm_yr']:+.2f} mm/yr, "
              f"NMAD {a['nmad_mm_yr']:.2f}, r = {a['pearson_r']:.2f}")

    if not rows_v:
        raise SystemExit(f"No matched products in {pts}: run the products step first")
    V, T, G, C = pd.DataFrame(rows_v), pd.DataFrame(rows_t), pd.DataFrame(rows_g), pd.DataFrame(cover)
    V.to_csv(out / "velocity_agreement.csv", index=False)
    T.to_csv(out / "timeseries_agreement.csv", index=False)
    G.to_csv(out / "elevation_gradient.csv", index=False)
    C.to_csv(out / "coverage.csv", index=False)
    figure_gradient(grad_fig, out / "fig_elevation_gradient.png", n_far, seed)
    scan = pd.read_csv(pts / "offset_scan.csv") if (pts / "offset_scan.csv").is_file() else pd.DataFrame()
    write_report(cfg, variant, out, C, V, T, G, scan, summary)
    print(f"  report: {out / 'comparison_report.md'}")
    return out


def _md(df: pd.DataFrame, floatfmt: str = "{:.2f}") -> str:
    if df is None or not len(df):
        return "_(none)_\n"
    cols = list(df.columns)
    lines = ["| " + " | ".join(cols) + " |", "|" + "|".join("---" for _ in cols) + "|"]
    for _, r in df.iterrows():
        cells = []
        for v in r:
            if isinstance(v, (float, np.floating)):
                cells.append("" if not np.isfinite(v) else floatfmt.format(v))
            else:
                cells.append("" if v is None else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines) + "\n"


def chain_text(cfg, variant: str) -> str:
    topo = cfg.mintpy.get("topographic_residual", True)
    if variant.startswith("hyp3_"):
        h = cfg.hyp3
        return (f"ASF HyP3 INSAR_ISCE_MULTI_BURST (ISCE2 run by ASF, {h.get('looks', '5x1')} looks; small-baseline "
                f"network of the {h.get('connections', 3)} nearest dates{' + one-year pairs' if h.get('annual_pairs', True) else ''}) "
                f"-> MintPy (troposphere: {variant[5:]}; DEM error correction: {topo}); UTM grid resampled "
                f"(nearest) to the SARscape lat/lon grid")
    return (f"ISCE2 topsStack (SLC) -> MiaplPy ({cfg.miaplpy.get('phase_linking', 'sequential_EMI')}, "
            f"{cfg.miaplpy.get('network', 'single_reference')} network) -> MintPy (troposphere: {variant}; "
            f"DEM error correction: {topo})")


def write_report(cfg, variant, out, C, V, T, G, scan, summary) -> None:
    best = pd.DataFrame()
    if len(scan):
        best = scan.loc[scan.groupby(["product", "positions"])["nmad"].idxmin()].copy()
        zero = scan[(scan.east_m == 0) & (scan.north_m == 0)][["product", "positions", "nmad"]]
        best = best.merge(zero.rename(columns={"nmad": "nmad_at_0"}), on=["product", "positions"], how="left")
        best = best[["product", "positions", "east_m", "north_m", "nmad", "nmad_at_0", "points"]]
    answers = cfg.export_dir / "answers" / f"answers_side_by_side_{variant}.md"
    vcols = ["product", "points", "median_difference_mm_yr", "nmad_mm_yr", "rmse_mm_yr", "pearson_r",
             "spearman_rho", "orthogonal_slope", "share_within_1_mm_yr", "share_within_2_mm_yr"]
    txt = [f"# SARscape E-SBAS vs open-source chain - variant '{variant}'",
           "",
           f"Created {datetime.now():%Y-%m-%d %H:%M}. Open-source chain: {chain_text(cfg, variant)}. "
           f"Common reference point of the open-source run: {cfg.mintpy.get('reference_lalo')}. "
           f"Both chains re-referenced to the **{cfg.compare.get('reference', 'headworks')}** median before comparison; "
           "velocities are least-squares rates over the common dates.",
           "", "## 1. Coverage", "", _md(C, "{:.3f}"),
           "## 2. Velocity agreement (open source - SARscape, same points, same dates)", "",
           _md(V[vcols]), "",
           "NMAD = 1.4826 x median absolute deviation of the difference (robust spread). "
           "orthogonal_slope = slope of the major axis of the open-source vs SARscape cloud (1 = same scale).",
           "", "![ascending](fig_ascending.png)", "![descending](fig_descending.png)",
           "![vertical](fig_vertical.png)", "![east](fig_east.png)", "",
           "## 3. Agreement of the displacement histories", "", _md(T),
           "The per-date median difference (common_mode_<product>.csv, lower-right panels) shows what differs "
           "across the whole window at once: residual atmosphere, reference handling or unwrapping.", "",
           "## 4. Elevation dependence of the velocity (far field)", "", _md(G),
           "![gradient](fig_elevation_gradient.png)", "",
           "This is the test of the manuscript's Section 4.2: if the height-dependent bias is residual "
           "tropospheric stratification of the SARscape run, it should shrink with the weather-model variant "
           "('pyaps' / 'gacos') of the open-source chain.", "",
           "## 5. Geolocation check (offset scan)", "", _md(best),
           "Shift of the sampling positions (m east / north) that minimises the spread of the velocity "
           "difference; nmad_at_0 is the spread without a shift. For the decomposed products the scan is made "
           "both at SARscape's lon/lat fields and at the shapefile geometry, which tests the ~22 m offset of the "
           "manuscript copies.", "",
           "## 6. Answers of the agent on both product sets", "",
           (f"See [{answers.name}](../../answers/{answers.name})." if answers.is_file() else
            "Not yet run: on Windows, run kgproc/answers.py with the agent's Python (README, step 8)."), ""]
    (out / "comparison_report.md").write_text("\n".join(txt), encoding="utf-8")
