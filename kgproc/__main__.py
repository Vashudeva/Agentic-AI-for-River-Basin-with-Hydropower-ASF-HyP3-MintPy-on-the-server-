"""Command line:  python -m kgproc [-c config.toml] <step> [--track ascending|descending] [--variant ...]

Route 'isce'  (Linux/WSL2, config.toml):        inventory aux isce miaplpy correct products compare
Route 'hyp3'  (Windows, config_hyp3.toml):      hyp3-search hyp3-submit hyp3-status hyp3-download hyp3-mintpy products compare
'all' runs the steps of the route given in the config file (route = "isce" | "hyp3").
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

from .config import TRACKS, load
from .util import StepError, free_gb, which

STEPS = ["check", "inventory", "aux", "isce", "miaplpy", "correct",
         "hyp3-search", "hyp3-submit", "hyp3-status", "hyp3-download", "hyp3-mintpy",
         "products", "compare", "all", "set-earthdata-token", "set-cds-key"]


def _trust_system_certificates() -> None:
    """Use the Windows/OS certificate store (needed behind university TLS inspection)."""
    try:
        import truststore
        truststore.inject_into_ssl()
    except Exception:                                   # noqa: BLE001
        pass


def cmd_check(cfg, _a) -> int:
    import importlib
    print(f"config    : {cfg.path}  (route: {cfg.route})")
    print(f"python    : {sys.executable} ({sys.version.split()[0]})")
    ok = True
    mods = (("mintpy", "hyp3_sdk", "asf_search", "h5py", "geopandas", "pyogrio", "shapely", "pyproj", "osgeo.gdal",
             "scipy", "pyaps3") if cfg.route == "hyp3" else
            ("isce", "mintpy", "miaplpy", "h5py", "geopandas", "pyogrio", "shapely", "pyproj", "osgeo.gdal", "scipy"))
    for mod in mods:
        try:
            m = importlib.import_module(mod)
            print(f"  {mod:<11}{getattr(m, '__version__', 'ok')}")
        except Exception as exc:                        # noqa: BLE001
            ok = False
            print(f"  {mod:<11}MISSING ({type(exc).__name__}: {exc})")
    if cfg.route != "hyp3":
        for exe in ("stackSentinel.py", "miaplpyApp", "smallbaselineApp.py", "geocode.py", "snaphu", "eof", "sardem"):
            w = which(exe)
            ok &= w is not None
            print(f"  {exe:<20}{w or 'NOT FOUND'}")
    print("\nInputs")
    for key, t in cfg.tracks.items():
        if cfg.route != "hyp3":
            print(f"  {key:<11} slc_dir {'ok' if t.slc_dir.is_dir() else 'NOT FOUND'}: {t.slc_dir}")
        print(f"  {key:<11} CG report {'ok' if t.sarscape_dates.is_file() else 'NOT FOUND'}: {t.sarscape_dates}")
    for key, p in cfg.sarscape_products.items():
        print(f"  SARscape {key:<11}{'ok' if p.is_file() else 'NOT FOUND'}: {p}")
    if cfg.agent_dir:
        s = cfg.agent_dir / "settings.toml"
        print(f"  agent settings {'ok' if s.is_file() else 'NOT FOUND'}: {s}")
    if cfg.route == "hyp3":
        from .credentials import load_token
        print(f"  Earthdata token {'ok' if load_token(cfg) else 'NOT SET - run set-earthdata-token'}")
    if "pyaps" in cfg.tropo:
        rc = Path.home() / ".cdsapirc"
        key_ok = rc.is_file()
        try:
            import configparser
            import pyaps3
            cp = configparser.RawConfigParser(delimiters="=")
            cp.read(Path(pyaps3.__file__).with_name("model.cfg"))
            key_ok |= cp.has_option("CDS", "key") and "your-personal" not in cp.get("CDS", "key")
        except Exception:                               # noqa: BLE001
            pass
        print(f"  ERA5 (CDS) key {'ok' if key_ok else 'NOT SET - run set-cds-key (needed for variant pyaps)'}")
    print(f"\nFree space: work_dir {free_gb(cfg.work_dir):.0f} GB, export_dir {free_gb(cfg.export_dir):.0f} GB")
    return 0 if ok else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(prog="python -m kgproc", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("step", choices=STEPS)
    ap.add_argument("-c", "--config", default=None, help="config file (default: config.toml next to kgproc)")
    ap.add_argument("--track", choices=TRACKS, help="only this stack (default: both)")
    ap.add_argument("--variant", help="only this variant (default: all in the config)")
    ap.add_argument("--allow-missing", action="store_true", help="inventory: continue if dates are missing")
    ap.add_argument("--force", action="store_true", help="isce: delete the ISCE2 folder and start again")
    ap.add_argument("--yes", action="store_true", help="hyp3-submit: really submit (spends HyP3 credits)")
    ap.add_argument("--wait", action="store_true", help="hyp3-download: wait for unfinished jobs first")
    a = ap.parse_args(argv)
    _trust_system_certificates()
    cfg = load(a.config)
    tracks = [a.track] if a.track else list(TRACKS)
    variants = [a.variant] if a.variant else cfg.variants
    step = a.step
    t0 = time.time()
    try:
        if step == "check":
            return cmd_check(cfg, a)
        if step == "set-earthdata-token":
            from .credentials import set_token
            set_token(cfg)
            return 0
        if step == "set-cds-key":
            from .credentials import set_cds_key
            set_cds_key()
            return 0

        hyp3_all = step == "all" and cfg.route == "hyp3"
        isce_all = step == "all" and cfg.route != "hyp3"
        # ---------------- route isce (Linux)
        if step == "inventory" or isce_all:
            from . import inventory
            for k in tracks:
                print(f"[inventory] {k}")
                inventory.run(cfg, k, allow_missing=a.allow_missing)
        if step == "aux" or isce_all:
            from . import auxdata
            for k in tracks:
                print(f"[aux] {k}")
                auxdata.orbits(cfg, k)
            auxdata.dem(cfg)
            auxdata.aux_dir(cfg)
        if step == "isce" or isce_all:
            from . import auxdata, isce_stack
            dem = auxdata.dem(cfg)
            for k in tracks:
                print(f"[isce] {k}")
                isce_stack.run_step(cfg, k, dem, force=a.force)
        if step == "miaplpy" or isce_all:
            from . import miaplpy_run
            for k in tracks:
                print(f"[miaplpy] {k}")
                miaplpy_run.run_miaplpy(cfg, k)
        if step == "correct" or isce_all:
            from . import miaplpy_run
            for k in tracks:
                for v in variants:
                    print(f"[correct] {k} / {v}")
                    miaplpy_run.run_correct(cfg, k, v)
        # ---------------- route hyp3 (Windows)
        if step == "hyp3-search" or hyp3_all:
            from . import hyp3_run
            for k in tracks:
                print(f"[hyp3-search] {k}")
                hyp3_run.search(cfg, k)
        if step == "hyp3-submit" or hyp3_all:
            from . import hyp3_run
            for k in tracks:
                print(f"[hyp3-submit] {k}")
                hyp3_run.submit(cfg, k, yes=a.yes)
            if hyp3_all and not a.yes:
                return 0
        if step == "hyp3-status":
            from . import hyp3_run
            for k in tracks:
                hyp3_run.status(cfg, k)
        if step == "hyp3-download" or hyp3_all:
            from . import hyp3_run
            for k in tracks:
                print(f"[hyp3-download] {k}")
                hyp3_run.download(cfg, k, wait=a.wait or hyp3_all)
        if step == "hyp3-mintpy" or hyp3_all:
            from . import hyp3_run
            tropo = [v.removeprefix("hyp3_") for v in variants]
            for k in tracks:
                for v in tropo:
                    print(f"[hyp3-mintpy] {k} / {v}")
                    hyp3_run.run_mintpy(cfg, k, v)
        # ---------------- shared
        if step == "products" or step == "all":
            from . import products
            for v in variants:
                print(f"[products] {v}")
                products.run_step(cfg, v)
        if step == "compare" or step == "all":
            from . import compare
            for v in variants:
                print(f"[compare] {v}")
                compare.run_step(cfg, v)
    except StepError as exc:
        print(f"\nERROR: {exc}")
        return 2
    print(f"\nfinished in {(time.time() - t0) / 3600:.2f} h")
    return 0


if __name__ == "__main__":
    sys.exit(main())
