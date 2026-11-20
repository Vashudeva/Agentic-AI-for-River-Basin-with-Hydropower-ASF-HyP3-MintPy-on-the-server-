"""Route 'hyp3': ASF HyP3 burst interferograms (ISCE2, run by ASF) + MintPy, on Windows.

Steps (python -m kgproc <step>, config_hyp3.toml):
  hyp3-search    Sentinel-1 bursts of the SARscape dates that cover the window (ASF catalogue)
  hyp3-submit    small-baseline network -> INSAR_ISCE_MULTI_BURST jobs (asks before spending credits)
  hyp3-status    progress of the jobs
  hyp3-download  download finished products, keep only the window (same grid for every pair)
  hyp3-mintpy    MintPy time series, once per tropospheric variant -> <work_dir>/<track>/geo/hyp3_<variant>
then the shared steps: products --variant hyp3_<v>, compare --variant hyp3_<v>, answers.py.

HyP3 runs ISCE2 on ASF's computers from ASF's copies of the same acquisitions; nothing is
processed from the local SAFE folders. Jobs submitted to HyP3 are publicly visible (ASF
policy), and finished products are deleted by ASF after about two weeks: download promptly.
"""
from __future__ import annotations

import json
import shutil
import sys
import zipfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
import pandas as pd

from .config import Config
from .inventory import parse_cg_report
from .util import StepError, run

LAYERS_PAIR = ("unw_phase", "corr", "conncomp")
LAYERS_GEOMETRY = ("dem", "lv_theta", "lv_phi", "water_mask")


# --------------------------------------------------------------------------- helpers
def job_name(cfg: Config, key: str) -> str:
    t = cfg.tracks[key]
    return f"{cfg.hyp3.get('job_name_prefix', 'kali')}_{key[:3]}_{t.relative_orbit}"


def window_wkt(cfg: Config, key: str) -> str:
    S, N, W, E = cfg.tracks[key].window_snwe
    return f"POLYGON(({W} {S},{E} {S},{E} {N},{W} {N},{W} {S}))"


def connect(cfg: Config):
    from hyp3_sdk import HyP3
    from .credentials import load_token
    tok = load_token(cfg)
    return HyP3(token=tok) if tok else HyP3(prompt="password")


def burst_date(granule: str) -> str:
    """S1_336254_IW2_20200211T123025_VV_8F2E-BURST -> 20200211"""
    return granule.split("_")[3][:8]


def _burst_number(full_id: str) -> int:
    return int(str(full_id).split("_")[1])


# --------------------------------------------------------------------------- search
def search(cfg: Config, key: str) -> pd.DataFrame:
    import asf_search as asf
    t = cfg.tracks[key]
    ref, dates, trk = parse_cg_report(t.sarscape_dates)
    if trk != t.relative_orbit:
        raise StepError(f"{t.sarscape_dates} is for track {trk}, config says {t.relative_orbit}")
    d0 = datetime.strptime(dates[0], "%Y%m%d") - timedelta(days=1)
    d1 = datetime.strptime(dates[-1], "%Y%m%d") + timedelta(days=2)
    print(f"  SARscape list: {len(dates)} dates {dates[0]} .. {dates[-1]}; searching ASF bursts "
          f"(track {t.relative_orbit}, IW{t.swath}, VV) over the window ...")
    res = asf.search(dataset=asf.DATASET.SLC_BURST, intersectsWith=window_wkt(cfg, key),
                     relativeOrbit=t.relative_orbit, polarization="VV", start=d0, end=d1)
    rows = []
    for r in res:
        p = r.properties
        b = p.get("burst", {})
        rows.append({"date": str(p["startTime"])[:10].replace("-", ""), "granule": p["sceneName"],
                     "full_burst_id": b.get("fullBurstID"), "subswath": b.get("subswath"),
                     "burst_index": b.get("burstIndex"), "start": p["startTime"]})
    df = pd.DataFrame(rows)
    if df.empty:
        raise StepError("ASF returned no bursts: check the window, track and internet access")
    df = df[(df.subswath == f"IW{t.swath}") & df.date.isin(dates)].drop_duplicates(["date", "full_burst_id"])

    # bursts of the window: those present on (nearly) every SARscape date
    share = df.groupby("full_burst_id").date.nunique() / len(dates)
    ids = sorted(share[share >= 0.9].index, key=_burst_number)
    nums = [_burst_number(i) for i in ids]
    if not ids:
        raise StepError("No burst covers the window on most dates")
    if nums != list(range(nums[0], nums[0] + len(nums))):
        raise StepError(f"Bursts are not contiguous: {ids}")
    df = df[df.full_burst_id.isin(ids)].sort_values(["date", "full_burst_id"])
    complete = df.groupby("date").full_burst_id.nunique() == len(ids)
    ok_dates = sorted(complete[complete].index)
    missing = sorted(set(dates) - set(ok_dates))
    out = cfg.hyp3_dir(key)
    out.mkdir(parents=True, exist_ok=True)
    df[df.date.isin(ok_dates)].to_csv(out / "bursts.csv", index=False)
    summary = {"track": key, "relative_orbit": t.relative_orbit, "swath": f"IW{t.swath}",
               "burst_ids": ids, "sarscape_dates": len(dates), "dates_with_all_bursts": len(ok_dates),
               "missing_dates": missing, "sarscape_reference": ref,
               "created": datetime.now().isoformat(timespec="seconds")}
    (out / "bursts_summary.json").write_text(json.dumps(summary, indent=1))
    print(f"  bursts: {', '.join(ids)}; {len(ok_dates)}/{len(dates)} dates complete"
          + (f"; missing: {', '.join(missing)}" if missing else ""))
    return df


def _bursts(cfg: Config, key: str) -> tuple[dict[str, list[str]], list[str]]:
    f = cfg.hyp3_dir(key) / "bursts.csv"
    if not f.is_file():
        raise StepError(f"Run hyp3-search first ({f} missing)")
    df = pd.read_csv(f, dtype=str)
    by_date = {d: sorted(g.granule, key=lambda s: int(s.split("_")[1])) for d, g in df.groupby("date")}
    return by_date, sorted(by_date)


# --------------------------------------------------------------------------- network
def network(dates: list[str], connections: int = 3, annual: bool = True,
            annual_days: tuple[int, int] = (330, 400)) -> list[tuple[str, str]]:
    """Each date with its `connections` next dates, plus for each date the pair closest to one
    year inside annual_days (if any)."""
    d = pd.to_datetime(dates, format="%Y%m%d")
    pairs = set()
    for i in range(len(dates)):
        for k in range(1, connections + 1):
            if i + k < len(dates):
                pairs.add((dates[i], dates[i + k]))
        if annual:
            dt = (d - d[i]).days.to_numpy()
            cand = [j for j in range(i + 1, len(dates)) if annual_days[0] <= dt[j] <= annual_days[1]]
            if cand:
                j = min(cand, key=lambda j: abs(dt[j] - 365))
                pairs.add((dates[i], dates[j]))
    return sorted(pairs)


def planned_pairs(cfg: Config, key: str) -> list[tuple[str, str]]:
    _, dates = _bursts(cfg, key)
    h = cfg.hyp3
    return network(dates, int(h.get("connections", 3)), bool(h.get("annual_pairs", True)),
                   tuple(h.get("annual_days", [330, 400])))


# --------------------------------------------------------------------------- jobs
def _existing(hyp3, name: str) -> dict[tuple[str, str], object]:
    out = {}
    for job in hyp3.find_jobs(name=name):
        p = job.job_parameters or {}
        refs, secs = p.get("reference") or [], p.get("secondary") or []
        if refs and secs:
            pair = (burst_date(refs[0]), burst_date(secs[0]))
            if pair not in out or out[pair].status_code == "FAILED":
                out[pair] = job
    return out


def _estimate(hyp3, prepared: list[dict]) -> float | None:
    """Credit cost via the API's validate-only submission (nothing is submitted)."""
    try:
        r = hyp3.session.post(hyp3._get_endpoint_url("/jobs"), json={"jobs": prepared[:200], "validate_only": True})
        r.raise_for_status()
        costs = [j.get("credit_cost") for j in r.json().get("jobs", [])]
        if costs and all(c is not None for c in costs):
            return float(np.mean(costs)) * len(prepared)
    except Exception:                                    # noqa: BLE001
        return None
    return None


def submit(cfg: Config, key: str, yes: bool = False) -> None:
    from hyp3_sdk import HyP3
    by_date, dates = _bursts(cfg, key)
    pairs = planned_pairs(cfg, key)
    name = job_name(cfg, key)
    hyp3 = connect(cfg)
    have = _existing(hyp3, name)
    todo = [p for p in pairs if p not in have or have[p].status_code == "FAILED"]
    looks = str(cfg.hyp3.get("looks", "5x1"))
    prepared = [HyP3.prepare_insar_isce_multi_burst_job(
        reference=by_date[a], secondary=by_date[b], name=name,
        apply_water_mask=bool(cfg.hyp3.get("apply_water_mask", False)), looks=looks) for a, b in todo]
    (cfg.hyp3_dir(key) / "pairs.csv").write_text("reference,secondary\n" + "\n".join(f"{a},{b}" for a, b in pairs) + "\n")
    print(f"  network: {len(pairs)} pairs over {len(dates)} dates ({cfg.hyp3.get('connections', 3)} nearest"
          f"{' + annual' if cfg.hyp3.get('annual_pairs', True) else ''}); already at HyP3: {len(pairs) - len(todo)}; "
          f"to submit: {len(todo)} (job name '{name}', looks {looks}, {len(by_date[dates[0]])} bursts each)")
    if not todo:
        return
    est = _estimate(hyp3, prepared)
    try:
        left = hyp3.check_credits()
    except Exception:                                    # noqa: BLE001
        left = None
    print(f"  estimated cost: {f'{est:.0f} credits' if est is not None else 'unknown (see hyp3.costs())'}; "
          f"credits left this month: {left if left is not None else 'unknown'}")
    if not yes:
        print("  Nothing submitted. Run again with --yes to submit these jobs.")
        return
    submitted = 0
    for k in range(0, len(prepared), 100):
        hyp3.submit_prepared_jobs(prepared[k:k + 100])
        submitted += len(prepared[k:k + 100])
        print(f"  submitted {submitted}/{len(prepared)}")
    print("  HyP3 now processes the jobs (typically minutes to hours each). Check with hyp3-status.")


def status(cfg: Config, key: str) -> dict:
    hyp3 = connect(cfg)
    jobs = list(hyp3.find_jobs(name=job_name(cfg, key)))
    counts = pd.Series([j.status_code for j in jobs]).value_counts().to_dict() if jobs else {}
    expired = sum(1 for j in jobs if j.succeeded() and j.expired())
    exp = [j.expiration_time for j in jobs if j.succeeded() and j.expiration_time]
    print(f"  {job_name(cfg, key)}: {len(jobs)} jobs {counts}"
          + (f"; {expired} expired" if expired else "")
          + (f"; first expiry {min(exp):%Y-%m-%d}" if exp else ""))
    return counts


# --------------------------------------------------------------------------- download + clip
def _grid(cfg: Config, key: str, sample_tif: Path) -> dict:
    """Projection, pixel size and bounds shared by every pair of the track."""
    f = cfg.hyp3_dir(key) / "grid.json"
    if f.is_file():
        return json.loads(f.read_text())
    from osgeo import gdal, osr
    from pyproj import Transformer
    ds = gdal.Open(str(sample_tif))
    srs = osr.SpatialReference(wkt=ds.GetProjection())
    srs.AutoIdentifyEPSG()
    epsg = int(srs.GetAuthorityCode(None))
    res = abs(ds.GetGeoTransform()[1])
    S, N, W, E = cfg.tracks[key].window_snwe
    lon = np.r_[np.linspace(W, E, 50), np.full(50, E), np.linspace(E, W, 50), np.full(50, W)]
    lat = np.r_[np.full(50, S), np.linspace(S, N, 50), np.full(50, N), np.linspace(N, S, 50)]
    x, y = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True).transform(lon, lat)
    pad = 5 * res
    b = [np.floor((min(x) - pad) / res) * res, np.floor((min(y) - pad) / res) * res,
         np.ceil((max(x) + pad) / res) * res, np.ceil((max(y) + pad) / res) * res]
    grid = {"epsg": epsg, "res": res, "bounds": [float(v) for v in b]}
    f.write_text(json.dumps(grid, indent=1))
    return grid


def _clip(src: Path, dst: Path, grid: dict) -> None:
    from osgeo import gdal
    gdal.UseExceptions()
    gdal.Warp(str(dst), str(src), dstSRS=f"EPSG:{grid['epsg']}", outputBounds=grid["bounds"],
              xRes=grid["res"], yRes=grid["res"], targetAlignedPixels=True, resampleAlg="near",
              creationOptions=["COMPRESS=DEFLATE", "TILED=YES"])


def _unpack(cfg: Config, key: str, zip_path: Path, need_geometry: bool) -> str:
    product = zip_path.stem
    out = cfg.hyp3_dir(key) / product
    tmp = cfg.hyp3_dir(key) / "_tmp" / product
    tmp.mkdir(parents=True, exist_ok=True)
    wanted = list(LAYERS_PAIR) + (list(LAYERS_GEOMETRY) if need_geometry else [])
    with zipfile.ZipFile(zip_path) as z:
        for m in z.namelist():
            base = Path(m).name
            if base == f"{product}.txt" or any(base == f"{product}_{lay}.tif" for lay in wanted):
                (tmp / base).write_bytes(z.read(m))
    if not (tmp / f"{product}_unw_phase.tif").is_file() or not (tmp / f"{product}.txt").is_file():
        raise StepError(f"{zip_path.name}: unwrapped phase or metadata missing")
    grid = _grid(cfg, key, tmp / f"{product}_unw_phase.tif")
    out.mkdir(parents=True, exist_ok=True)
    for tif in tmp.glob("*.tif"):
        _clip(tif, out / f"{tif.stem}_clip.tif", grid)
    shutil.copy2(tmp / f"{product}.txt", out / f"{product}.txt")
    shutil.rmtree(tmp, ignore_errors=True)
    return product


def download(cfg: Config, key: str, wait: bool = False) -> None:
    name = job_name(cfg, key)
    hyp3 = connect(cfg)
    batch = hyp3.find_jobs(name=name)
    if wait:
        pending = [j for j in batch if not j.complete()]
        if pending:
            print(f"  waiting for {len(pending)} unfinished job(s) ...")
            from hyp3_sdk import Batch
            hyp3.watch(Batch(pending), timeout=int(cfg.hyp3.get("wait_hours", 24)) * 3600, interval=300)
            batch = hyp3.find_jobs(name=name)
    root = cfg.hyp3_dir(key)
    done_file = root / "downloaded.json"
    done = json.loads(done_file.read_text()) if done_file.is_file() else {}
    have_geometry = any(root.glob("*/*_dem_clip.tif"))
    n_new = n_exp = 0
    for job in batch:
        if not job.succeeded():
            continue
        if job.job_id in done:
            continue
        if job.expired():
            n_exp += 1
            continue
        zips = job.download_files(root / "_zips")
        product = _unpack(cfg, key, Path(zips[0]), need_geometry=not have_geometry)
        have_geometry = True
        if not cfg.hyp3.get("keep_zips", False):
            for z in zips:
                Path(z).unlink(missing_ok=True)
        done[job.job_id] = product
        done_file.write_text(json.dumps(done, indent=1))
        n_new += 1
        print(f"    {product}")
    counts = pd.Series([j.status_code for j in batch]).value_counts().to_dict() if len(batch) else {}
    print(f"  {name}: {n_new} new product(s), {len(done)} in total; jobs {counts}"
          + (f"; {n_exp} expired before download (resubmit them with hyp3-submit --yes)" if n_exp else ""))


# --------------------------------------------------------------------------- MintPy
def mintpy_template(cfg: Config, key: str, variant: str) -> Path:
    from .miaplpy_run import _kv
    h = cfg.hyp3_dir(key).as_posix()
    m, e = cfg.mintpy, cfg.export
    lat, lon = m["reference_lalo"]
    opts = [
        ("mintpy.load.processor", "hyp3"),
        ("mintpy.load.unwFile", f"{h}/*/*_unw_phase_clip.tif"),
        ("mintpy.load.corFile", f"{h}/*/*_corr_clip.tif"),
        ("mintpy.load.connCompFile", f"{h}/*/*_conncomp_clip.tif"),
        ("mintpy.load.demFile", f"{h}/*/*_dem_clip.tif"),
        ("mintpy.load.incAngleFile", f"{h}/*/*_lv_theta_clip.tif"),
        ("mintpy.load.azAngleFile", f"{h}/*/*_lv_phi_clip.tif"),
        ("mintpy.load.waterMaskFile", f"{h}/*/*_water_mask_clip.tif"),
        ("mintpy.compute.maxMemory", m.get("max_memory_gb", 8)),
        ("mintpy.reference.lalo", f"{lat}, {lon}"),
        ("mintpy.network.coherenceBased", "no"),
        ("mintpy.unwrapError.method", cfg.hyp3.get("unwrap_error_correction", "no")),
        ("mintpy.networkInversion.minTempCoh", e.get("min_tcoh", 0.7)),
        ("mintpy.troposphericDelay.method", variant),
        ("mintpy.topographicResidual", "yes" if m.get("topographic_residual", True) else "no"),
        ("mintpy.deramp", m.get("deramp", "no")),
        ("mintpy.geocode", "no"),                         # HyP3 products are already geocoded (UTM)
        ("mintpy.save.kmz", "no"),
        ("mintpy.save.hdfEos5", "no"),
        ("mintpy.plot", "no"),
    ]
    if variant == "pyaps":
        wdir = m.get("weather_dir") or (cfg.work_dir / "weather").as_posix()
        Path(wdir).mkdir(parents=True, exist_ok=True)
        opts += [("mintpy.troposphericDelay.weatherModel", "ERA5"), ("mintpy.troposphericDelay.weatherDir", wdir)]
    elif variant == "gacos":
        if not m.get("gacos_dir"):
            raise StepError("troposphere variant 'gacos' needs mintpy.gacos_dir")
        opts.append(("mintpy.troposphericDelay.gacosDir", str(m["gacos_dir"]).replace("\\", "/")))
    path = cfg.track_dir(key) / f"mintpy_hyp3_{key}_{variant}.cfg"
    path.write_text(f"# MintPy options, HyP3 route, variant '{variant}' (written by kgproc)\n" + _kv(opts))
    return path


def run_mintpy(cfg: Config, key: str, variant: str) -> Path:
    from .util import done, mark
    work = cfg.hyp3_mintpy_dir(key)
    out = cfg.geo_dir(key, f"hyp3_{variant}")
    marker = cfg.marker(key, f"hyp3_mintpy.{variant}")
    if done(marker) and out.is_dir():
        return out
    pairs = list(cfg.hyp3_dir(key).glob("*/*_unw_phase_clip.tif"))
    if len(pairs) < 10:
        raise StepError(f"Only {len(pairs)} interferograms in {cfg.hyp3_dir(key)}: run hyp3-download first")
    tmpl = mintpy_template(cfg, key, variant)
    work.mkdir(parents=True, exist_ok=True)
    log = cfg.track_dir(key) / "logs" / f"mintpy_hyp3_{variant}.log"
    cmd = [sys.executable, "-m", "mintpy.cli.smallbaselineApp", str(tmpl), "--dir", str(work)]
    if (work / "timeseries.h5").is_file():
        cmd += ["--start", "correct_LOD"]                # inversion done: only corrections for this variant
    print(f"  MintPy {key}/{variant}: {len(pairs)} interferograms (log: {log})")
    run(cmd, log, cwd=work)

    from mintpy.smallbaselineApp import TimeSeriesAnalysis
    from mintpy.utils import readfile, utils as ut
    template = ut.check_template_auto_value(readfile.read_template(str(work / "smallbaselineApp.cfg")))
    template["mintpy.troposphericDelay.method"] = False if variant == "no" else variant
    ts = Path(TimeSeriesAnalysis.get_timeseries_filename(template, str(work))["velocity"]["input"])
    if not ts.is_absolute():
        ts = work / ts
    if not ts.is_file():
        raise StepError(f"Expected time series {ts} not found; see {log}")
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ts, out / "geo_timeseries.h5")
    for src, dst in (("temporalCoherence.h5", "geo_temporalCoherence.h5"), ("velocity.h5", "geo_velocity.h5"),
                     ("inputs/geometryGeo.h5", "geo_geometryRadar.h5")):
        if (work / src).is_file():
            shutil.copy2(work / src, out / dst)
    (out / "SOURCE.txt").write_text(f"HyP3 route\ntime series: {ts.name}\ntemplate: {tmpl.name}\n")
    print(f"  {key}/hyp3_{variant}: {ts.name} -> {out}")
    mark(marker, ts.name)
    return out
