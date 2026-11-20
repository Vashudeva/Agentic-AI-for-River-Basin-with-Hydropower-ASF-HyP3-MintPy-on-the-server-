"""Step 'inventory': find the Sentinel-1 files of the SARscape dates and link them per track.

The acquisition list is read from SARscape's connection-graph report (CG_report.txt), so
the open-source run uses exactly the images of the E-SBAS run. Each SAFE file is checked
for relative orbit, pass direction, VV polarisation and coverage of the analysis window.
The selected files are linked into <work_dir>/<track>/slc, which is what ISCE2 reads.
"""
from __future__ import annotations

import json
import re
import zipfile
from datetime import datetime
from pathlib import Path

from shapely.geometry import Polygon, box
from shapely.ops import unary_union

from .config import Config, Track

NAME = re.compile(r"^(S1[ABCD])_IW_SLC__1S([SD])(V|H)_(\d{8})T(\d{6})_(\d{8})T(\d{6})_(\d{6})_")
CG_LINE = re.compile(r"^(REFERENCE|SECONDARY):\s*\d+\s*\[sentinel1_(\d+)_(\d{8})_")


def parse_cg_report(path: Path) -> tuple[str, list[str], int]:
    """Reference date, all dates (sorted) and track number from SARscape's CG_report.txt."""
    ref, dates, track = None, [], None
    for line in path.read_text(errors="replace").splitlines():
        m = CG_LINE.match(line.strip())
        if m:
            kind, trk, date = m.groups()
            track = int(trk)
            dates.append(date)
            if kind == "REFERENCE":
                ref = date
    if not dates:
        raise ValueError(f"No acquisitions found in {path}")
    return ref, sorted(set(dates)), track


def _manifest_text(path: Path) -> str:
    if path.suffix.lower() == ".zip":
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.endswith("manifest.safe")]
            if not names:
                raise ValueError("manifest.safe not found")
            return z.read(sorted(names, key=len)[0]).decode("utf-8", errors="replace")
    return (path / "manifest.safe").read_text(errors="replace")


def read_safe(path: Path) -> dict:
    """Metadata of one Sentinel-1 IW SLC product (.zip or .SAFE)."""
    m = NAME.match(path.name)
    if not m:
        raise ValueError("not a Sentinel-1 IW SLC name")
    mission, _, pol, d0, t0, _, _, abs_orbit = m.groups()
    info = {"file": str(path), "name": path.name, "mission": mission,
            "vv": pol == "V", "date": d0, "start": f"{d0}T{t0}", "absolute_orbit": int(abs_orbit)}
    text = _manifest_text(path)
    ro = re.search(r'relativeOrbitNumber type="start">(\d+)<', text)
    ps = re.search(r"<s1:pass>(\w+)</s1:pass>", text)
    fp = re.search(r"<gml:coordinates>([^<]+)</gml:coordinates>", text)
    info["relative_orbit"] = int(ro.group(1)) if ro else None
    info["pass"] = ps.group(1).upper() if ps else None
    if fp:
        pts = [tuple(float(v) for v in pair.split(",")) for pair in fp.group(1).split()]
        info["footprint"] = [(lon, lat) for lat, lon in pts]          # manifest gives lat,lon
    else:
        info["footprint"] = None
    return info


def scan(track: Track) -> list[dict]:
    pattern = "S1*_IW_SLC__1S*"
    found = track.slc_dir.rglob(pattern) if track.slc_recursive else track.slc_dir.glob(pattern)
    out = []
    for p in sorted(found):
        if not (p.suffix.lower() == ".zip" or (p.suffix.upper() == ".SAFE" and p.is_dir())):
            continue
        try:
            out.append(read_safe(p))
        except (ValueError, zipfile.BadZipFile, OSError) as exc:
            out.append({"file": str(p), "name": p.name, "error": str(exc)})
    return out


def run(cfg: Config, key: str, allow_missing: bool = False) -> dict:
    track = cfg.tracks[key]
    out_dir = cfg.track_dir(key) / "inventory"
    out_dir.mkdir(parents=True, exist_ok=True)
    if not track.slc_dir.is_dir():
        raise FileNotFoundError(f"slc_dir not found: {track.slc_dir}")

    ref, wanted, trk = parse_cg_report(track.sarscape_dates)
    if trk != track.relative_orbit:
        raise ValueError(f"{track.sarscape_dates} is for track {trk}, config says {track.relative_orbit}")
    if ref != track.reference_date:
        print(f"  note: SARscape reference date is {ref}, config reference_date is {track.reference_date}")
    print(f"  SARscape list: {len(wanted)} dates, {wanted[0]} .. {wanted[-1]}, reference {ref}")

    S, N, W, E = track.window_snwe
    window = box(W, S, E, N)
    files = scan(track)
    bad = [f for f in files if "error" in f]
    ours = [f for f in files if "error" not in f and f["relative_orbit"] == track.relative_orbit
            and f["pass"] == track.pass_direction and f["vv"]]
    by_date: dict[str, list[dict]] = {}
    for f in ours:
        by_date.setdefault(f["date"], []).append(f)

    rows, missing, partial, extra = [], [], [], sorted(set(by_date) - set(wanted))
    for d in wanted:
        fs = by_date.get(d, [])
        if not fs:
            missing.append(d)
            rows.append({"date": d, "status": "missing", "files": []})
            continue
        polys = [Polygon(f["footprint"]) for f in fs if f["footprint"]]
        covered = bool(polys) and unary_union(polys).buffer(1e-6).contains(window)
        use = [f for f in fs if not f["footprint"] or Polygon(f["footprint"]).intersects(window)]
        if not covered:
            partial.append(d)
        rows.append({"date": d, "status": "ok" if covered else "partial",
                     "mission": sorted({f["mission"] for f in use}),
                     "files": [f["file"] for f in use]})

    kinds = {Path(f).suffix.lower() for r in rows for f in r["files"]}
    if len(kinds) > 1:
        raise ValueError("Mixed .zip and .SAFE inputs: ISCE2 reads only one kind per folder. "
                         "Keep all scenes as .zip (or all as .SAFE).")

    summary = {"track": key, "relative_orbit": track.relative_orbit, "pass": track.pass_direction,
               "sarscape_reference": ref, "dates_expected": len(wanted),
               "dates_found": len(wanted) - len(missing), "missing": missing,
               "partial_coverage": partial, "other_dates_in_folder": extra,
               "unreadable_files": [b["file"] for b in bad],
               "missions": sorted({m for r in rows for m in r.get("mission", [])}),
               "created": datetime.now().isoformat(timespec="seconds"), "dates": rows}
    (out_dir / "inventory.json").write_text(json.dumps(summary, indent=1))
    (out_dir / "dates.txt").write_text("\n".join(wanted) + "\n")

    print(f"  found {summary['dates_found']}/{len(wanted)} dates"
          + (f"; missing {', '.join(missing)}" if missing else "")
          + (f"; window not fully covered on {', '.join(partial)}" if partial else "")
          + (f"; {len(bad)} unreadable file(s)" if bad else ""))
    if extra:
        print(f"  {len(extra)} other date(s) of this track in the folder are ignored (not in the SARscape list)")
    if (missing or partial) and not allow_missing:
        raise SystemExit("Inventory incomplete. Add the missing scenes (or run with --allow-missing "
                         "to process the dates that are available; the comparison then uses common dates).")

    link_dir = cfg.track_dir(key) / "slc"
    link_dir.mkdir(parents=True, exist_ok=True)
    for old in link_dir.iterdir():
        if old.is_symlink():
            old.unlink()
    n = 0
    for r in rows:
        if r["status"] == "missing":
            continue
        for f in r["files"]:
            (link_dir / Path(f).name).symlink_to(f)
            n += 1
    print(f"  linked {n} file(s) into {link_dir}")
    return summary
