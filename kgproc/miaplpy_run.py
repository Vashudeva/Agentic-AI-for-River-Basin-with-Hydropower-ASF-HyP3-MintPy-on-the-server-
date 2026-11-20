"""Steps 'miaplpy' and 'correct'.

miaplpy : MiaplPy from loading the coregistered SLCs to the network inversion
          (phase linking of PS and DS at full resolution, single-reference interferograms
          on the SARscape reference date, SNAPHU, MintPy network inversion).
correct : MintPy corrections and geocoding, once per tropospheric variant
          (mintpy.troposphere in config.toml). Geocoded files of each variant are copied
          to <work_dir>/<track>/geo/<variant>/ so that variants do not overwrite each other.
"""
from __future__ import annotations

import shutil
import time
from pathlib import Path

from .config import Config
from .util import StepError, done, mark, run, which

TROPO_TOKEN = {"no": None, "pyaps": "ERA5", "gacos": "GACOS", "height_correlation": "tropHgt"}


def _kv(lines: list[tuple[str, object]]) -> str:
    width = max(len(k) for k, _ in lines)
    return "\n".join(f"{k:<{width}} = {v}" for k, v in lines) + "\n"


def mintpy_options(cfg: Config, variant: str) -> list[tuple[str, object]]:
    m, e = cfg.mintpy, cfg.export
    lat, lon = m["reference_lalo"]
    S, N, W, E = m["geocode_snwe"]
    step = float(m.get("geocode_step_deg", 0.000125))
    opts = [
        ("mintpy.compute.maxMemory", m.get("max_memory_gb", 8)),
        ("mintpy.reference.lalo", f"{lat}, {lon}"),
        ("mintpy.network.coherenceBased", "no"),        # keep every single-reference interferogram
        ("mintpy.unwrapError.method", "no"),
        ("mintpy.networkInversion.minTempCoh", e.get("min_tcoh", 0.7)),
        ("mintpy.troposphericDelay.method", variant),
        ("mintpy.topographicResidual", "yes" if m.get("topographic_residual", True) else "no"),
        ("mintpy.deramp", m.get("deramp", "no")),
        ("mintpy.geocode", "yes"),
        ("mintpy.geocode.SNWE", f"{S}, {N}, {W}, {E}"),
        ("mintpy.geocode.laloStep", f"{-step}, {step}"),
        ("mintpy.geocode.interpMethod", "nearest"),
        ("mintpy.save.kmz", "no"),
        ("mintpy.save.hdfEos5", "no"),
        ("mintpy.plot", "no"),
    ]
    if variant == "pyaps":
        wdir = m.get("weather_dir") or str(cfg.work_dir / "weather")
        Path(wdir).mkdir(parents=True, exist_ok=True)
        opts += [("mintpy.troposphericDelay.weatherModel", "ERA5"),
                 ("mintpy.troposphericDelay.weatherDir", wdir)]
    elif variant == "gacos":
        if not m.get("gacos_dir"):
            raise StepError("troposphere variant 'gacos' needs mintpy.gacos_dir")
        opts.append(("mintpy.troposphericDelay.gacosDir", m["gacos_dir"]))
    elif variant not in TROPO_TOKEN:
        raise StepError(f"Unknown troposphere variant '{variant}'. Use one of {list(TROPO_TOKEN)}")
    return opts


def write_template(cfg: Config, key: str) -> Path:
    t, mp = cfg.tracks[key], cfg.miaplpy
    isce, work = cfg.isce_dir(key), cfg.miaplpy_dir(key)
    geom = isce / "merged" / "geom_reference"
    net = mp.get("network", "single_reference")
    ifg = work / "inverted" / f"interferograms_{net}"
    S, N, W, E = t.window_snwe
    lines = [
        ("miaplpy.load.processor", "isce"),
        ("miaplpy.load.updateMode", "auto"),
        ("miaplpy.load.compression", "auto"),
        ("miaplpy.load.autoPath", "no"),
        ("miaplpy.load.slcFile", f"{isce}/merged/SLC/*/*.slc.full"),
        ("miaplpy.load.metaFile", f"{isce}/reference/IW*.xml"),
        ("miaplpy.load.baselineDir", f"{isce}/baselines"),
        ("miaplpy.load.demFile", f"{geom}/hgt.rdr.full"),
        ("miaplpy.load.lookupYFile", f"{geom}/lat.rdr.full"),
        ("miaplpy.load.lookupXFile", f"{geom}/lon.rdr.full"),
        ("miaplpy.load.incAngleFile", f"{geom}/los.rdr.full"),
        ("miaplpy.load.azAngleFile", f"{geom}/los.rdr.full"),
        ("miaplpy.load.shadowMaskFile", f"{geom}/shadowMask.rdr.full" if (geom / "shadowMask.rdr.full").exists() else "None"),
        ("miaplpy.load.waterMaskFile", "None"),
        ("miaplpy.load.unwFile", f"{ifg}/*/filt_fine.unw"),
        ("miaplpy.load.corFile", f"{ifg}/*/filt_fine.cor"),
        ("miaplpy.load.connCompFile", f"{ifg}/*/filt_fine.unw.conncomp"),
        ("miaplpy.subset.lalo", f"{S}:{N},{W}:{E}"),
        ("miaplpy.multiprocessing.numProcessor", mp.get("num_processor", 4)),
        ("miaplpy.inversion.patchSize", mp.get("patch_size", 200)),
        ("miaplpy.inversion.ministackSize", mp.get("ministack_size", 10)),
        ("miaplpy.inversion.rangeWindow", mp.get("range_window", 19)),
        ("miaplpy.inversion.azimuthWindow", mp.get("azimuth_window", 5)),
        ("miaplpy.inversion.shpTest", mp.get("shp_test", "ks")),
        ("miaplpy.inversion.phaseLinkingMethod", mp.get("phase_linking", "sequential_EMI")),
        ("miaplpy.interferograms.networkType", net),
        ("miaplpy.interferograms.referenceDate", t.reference_date),
        ("miaplpy.timeseries.minTempCoh", mp.get("min_temp_coh", 0.5)),
    ] + mintpy_options(cfg, "no")
    work.mkdir(parents=True, exist_ok=True)
    path = cfg.track_dir(key) / f"kali_{key}.txt"
    path.write_text(f"# MiaplPy + MintPy options for the {key} stack (written by kgproc; edit config.toml instead)\n"
                    + _kv(lines))
    return path


def run_miaplpy(cfg: Config, key: str) -> None:
    if which("miaplpyApp") is None or which("snaphu") is None:
        raise StepError("miaplpyApp or snaphu not on PATH: run setup_env.sh")
    template = write_template(cfg, key)
    log = cfg.track_dir(key) / "logs" / "miaplpy.log"
    m = cfg.marker(key, "miaplpy")
    if done(m):
        print("  MiaplPy already finished (delete the marker in _done to rerun)")
        return
    run(["miaplpyApp", str(template), "--dir", str(cfg.miaplpy_dir(key)), "--end", "invert_network"], log,
        cwd=cfg.track_dir(key))
    ts = cfg.network_dir(key) / "timeseries.h5"
    if not ts.is_file():
        raise StepError(f"MiaplPy finished without {ts}; see {log}")
    mark(m)


def run_correct(cfg: Config, key: str, variant: str) -> Path:
    """MintPy corrections + geocoding for one variant; returns the folder with geocoded files."""
    net, out = cfg.network_dir(key), cfg.geo_dir(key, variant)
    m = cfg.marker(key, f"correct.{variant}")
    if done(m) and out.is_dir():
        return out
    tmpl = cfg.track_dir(key) / f"mintpy_{key}_{variant}.cfg"
    tmpl.write_text(f"# MintPy options, variant '{variant}' (written by kgproc)\n" + _kv(mintpy_options(cfg, variant)))
    log = cfg.track_dir(key) / "logs" / f"mintpy_{variant}.log"
    t0 = time.time()
    run(["smallbaselineApp.py", str(tmpl), "--start", "correct_LOD", "--dir", str(net)], log, cwd=net)

    geo = net / "geo"
    token = TROPO_TOKEN[variant]
    cands = [p for p in geo.glob("geo_timeseries*.h5") if p.stat().st_mtime >= t0 - 1]
    cands = [p for p in cands if (token in p.name if token else not any(
        tk in p.name for tk in TROPO_TOKEN.values() if tk))]
    if not cands:
        raise StepError(f"No geocoded time series for variant '{variant}' in {geo}; see {log}")
    ts = max(cands, key=lambda p: len(p.name))
    out.mkdir(parents=True, exist_ok=True)
    shutil.copy2(ts, out / "geo_timeseries.h5")
    for name in ("geo_temporalCoherence.h5", "geo_geometryRadar.h5", "geo_velocity.h5", "geo_maskTempCoh.h5"):
        if (geo / name).is_file():
            shutil.copy2(geo / name, out / name)
    (out / "SOURCE.txt").write_text(f"time series: {ts.name}\ntemplate: {tmpl.name}\n")

    ps = cfg.miaplpy_dir(key) / "maskPS.h5"
    lookup = net / "inputs" / "geometryRadar.h5"
    if ps.is_file() and lookup.is_file():
        run(["geocode.py", str(ps), "-l", str(lookup), "-t", str(tmpl), "--outdir", str(out)], log, cwd=net)
    print(f"  {key}/{variant}: {ts.name} -> {out}")
    mark(m, ts.name)
    return out
