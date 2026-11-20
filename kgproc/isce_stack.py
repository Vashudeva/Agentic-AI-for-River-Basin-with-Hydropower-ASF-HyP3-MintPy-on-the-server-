"""Step 'isce': ISCE2 topsStack, SLC workflow, restricted to the bursts of the window.

stackSentinel.py writes numbered run files; each is executed in order, its lines in
parallel (isce.num_proc). A finished run file is recorded, so the step resumes where it
stopped. Output used later: merged/SLC/<date>/<date>.slc.full, merged/geom_reference/*,
reference/IW*.xml and baselines/.
"""
from __future__ import annotations

import shutil
from pathlib import Path

from .auxdata import aux_dir
from .config import Config
from .util import StepError, done, mark, run, run_file, which


def command(cfg: Config, key: str, dem: Path) -> list[str]:
    t = cfg.tracks[key]
    S, N, W, E = t.window_snwe
    c = cfg.isce
    return ["stackSentinel.py",
            "-s", str(cfg.track_dir(key) / "slc"),
            "-o", str(cfg.track_dir(key) / "orbits"),
            "-a", str(aux_dir(cfg)),
            "-d", str(dem),
            "-w", str(cfg.isce_dir(key)),
            "-p", str(c.get("polarization", "vv")),
            "-W", "slc",
            "-n", str(t.swath),
            "-b", f"{S} {N} {W} {E}",
            "-C", str(c.get("coregistration", "NESD")),
            "-m", t.reference_date,
            "--num_proc", str(c.get("num_proc", 4)),
            "--num_proc4topo", str(c.get("num_proc", 4))]


def run_step(cfg: Config, key: str, dem: Path, force: bool = False) -> None:
    work = cfg.isce_dir(key)
    if which("stackSentinel.py") is None:
        raise StepError("stackSentinel.py not found: is ISCE2 installed in this environment?")
    log_dir = cfg.track_dir(key) / "logs" / "isce"
    if force and work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    runs = work / "run_files"
    if not runs.is_dir() or not any(runs.glob("run_*")):
        run(command(cfg, key, dem), log_dir / "stackSentinel.log", cwd=work)
    files = sorted(p for p in runs.glob("run_*") if not p.suffix and p.is_file())
    if not files:
        raise StepError(f"stackSentinel.py wrote no run files in {runs}; see {log_dir}")
    workers = int(cfg.isce.get("num_proc", 4))
    for f in files:
        m = cfg.marker(key, f"isce.{f.name}")
        if done(m):
            continue
        run_file(f, log_dir, workers, cwd=work)
        mark(m)
    check(cfg, key)
    if cfg.isce.get("cleanup"):
        for sub in ("coreg_secondarys", "secondarys", "ESD", "misreg"):
            shutil.rmtree(work / sub, ignore_errors=True)


def check(cfg: Config, key: str) -> None:
    import json
    work = cfg.isce_dir(key)
    inv = json.loads((cfg.track_dir(key) / "inventory" / "inventory.json").read_text())
    want = [r["date"] for r in inv["dates"] if r["status"] != "missing"]
    have = sorted(p.parent.name for p in (work / "merged" / "SLC").glob("*/*.slc.full"))
    geom = [work / "merged" / "geom_reference" / f"{n}.rdr.full" for n in ("hgt", "lat", "lon", "los")]
    lacking = sorted(set(want) - set(have))
    if lacking or not all(g.is_file() for g in geom):
        raise StepError(f"ISCE2 stack incomplete in {work}: missing SLC dates {lacking[:5]}"
                        f"{'...' if len(lacking) > 5 else ''}, geometry ok = {all(g.is_file() for g in geom)}")
    print(f"  ISCE2 stack ok: {len(have)} coregistered SLCs, geometry in merged/geom_reference")
