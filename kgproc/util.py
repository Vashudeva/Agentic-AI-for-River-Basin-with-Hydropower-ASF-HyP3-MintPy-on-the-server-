"""Shell helpers: environment for ISCE2/MiaplPy/MintPy, logged commands, run files."""
from __future__ import annotations

import os
import shlex
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path


class StepError(RuntimeError):
    pass


def tool_env() -> dict[str, str]:
    """Environment in which ISCE2 topsStack, MiaplPy and MintPy scripts run.

    Works without 'activating' the conda environment: everything is derived from the
    Python that runs kgproc (sys.prefix).
    """
    prefix = Path(sys.prefix)
    site = Path(next(p for p in sys.path if p.endswith("site-packages")))
    env = dict(os.environ)
    paths = [prefix / "bin", prefix / "share" / "isce2" / "topsStack", site / "isce" / "applications"]
    if os.name == "nt":                          # conda layout on Windows
        paths = [prefix, prefix / "Scripts", prefix / "Library" / "bin"] + paths
    env["PATH"] = os.pathsep.join([str(p) for p in paths] + [env.get("PATH", "")])
    py = [prefix / "share" / "isce2",            # 'topsStack' package used by the run files
          site / "isce" / "components",          # isceobj, contrib (MiaplPy unwrap_ifgram)
          site / "miaplpy"]                      # MiaplPy scripts import their siblings directly
    old = env.get("PYTHONPATH")
    env["PYTHONPATH"] = os.pathsep.join([str(p) for p in py] + ([old] if old else []))
    env.setdefault("ISCE_HOME", str(site / "isce"))
    share = prefix / "Library" / "share" if os.name == "nt" else prefix / "share"
    env.setdefault("PROJ_DATA", str(share / "proj"))
    env.setdefault("PROJ_LIB", str(share / "proj"))
    env.setdefault("GDAL_DATA", str(share / "gdal"))
    env.setdefault("OMP_NUM_THREADS", "1")       # parallelism comes from the run files / workers
    env["PYTHONUNBUFFERED"] = "1"
    return env


def which(name: str) -> str | None:
    return shutil.which(name, path=tool_env()["PATH"])


def run(cmd: list[str] | str, log: Path, cwd: Path | None = None, env: dict | None = None,
        check: bool = True) -> int:
    """Run a command, append stdout/stderr to a log file, echo a short status line."""
    log.parent.mkdir(parents=True, exist_ok=True)
    text = cmd if isinstance(cmd, str) else " ".join(shlex.quote(str(c)) for c in cmd)
    t0 = time.time()
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(f"\n$ {text}\n")
        fh.flush()
        proc = subprocess.run(text if isinstance(cmd, str) else [str(c) for c in cmd],
                              shell=isinstance(cmd, str), cwd=cwd, env=env or tool_env(),
                              stdout=fh, stderr=subprocess.STDOUT)
    dt = time.time() - t0
    status = "ok" if proc.returncode == 0 else f"FAILED ({proc.returncode})"
    print(f"    {status:<11}{dt/60:7.1f} min  {text[:110]}", flush=True)
    if check and proc.returncode != 0:
        raise StepError(f"Command failed, see {log}:\n  {text}")
    return proc.returncode


def run_file(path: Path, log_dir: Path, workers: int, cwd: Path) -> None:
    """Execute an ISCE2/MiaplPy run file: one command per line, lines independent."""
    lines = [ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]
    print(f"  {path.name}: {len(lines)} command(s)", flush=True)
    log = log_dir / f"{path.name}.log"

    def one(item):
        i, line = item
        return run(line, log_dir / f"{path.name}.{i:04d}.log", cwd=cwd, check=False), line

    with ThreadPoolExecutor(max_workers=max(1, workers)) as pool:
        results = list(pool.map(one, enumerate(lines)))
    failed = [ln for rc, ln in results if rc != 0]
    with open(log, "a", encoding="utf-8") as fh:
        fh.write(f"{len(lines) - len(failed)}/{len(lines)} commands succeeded\n")
    if failed:
        raise StepError(f"{len(failed)} command(s) of {path.name} failed; logs in {log_dir}\n  first: {failed[0]}")


def done(marker: Path) -> bool:
    return marker.is_file()


def mark(marker: Path, text: str = "") -> None:
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(time.strftime("%Y-%m-%d %H:%M:%S ") + text + "\n")


def free_gb(path: Path) -> float:
    p = path
    while not p.exists():
        p = p.parent
    return shutil.disk_usage(p).free / 1e9
