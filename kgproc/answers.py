"""The agent's answers on SARscape and on open-source products, side by side.

Run with the AGENT's Python on Windows (it needs only the agent and its packages):

    G:\\InSAR\\Kaligandaki\\Langgraph\\.venv\\Scripts\\python.exe kgproc\\answers.py ^
        --agent-dir G:\\InSAR\\Kaligandaki\\Langgraph\\insar-hydro-agent\\agent ^
        --settings-dir G:\\InSAR\\Kaligandaki\\Langgraph\\opensource\\agent_settings ^
        --variant no

Every settings file (sarscape, sarscape_matched_<variant>, opensource_at_sarscape_<variant>,
opensource_native_<variant>) is run in its own process with KALI_SETTINGS pointing to it,
through the agent's own functions (screening.build, elevation.summary, sites.zone_table),
i.e. exactly the numbers the agent and the manuscript report. Each run writes into its own
work_root (<export_dir>/agent_runs/<name>), so the live agent cache is not touched.

Output: <export_dir>/answers/answers_<run>.json and answers_side_by_side_<variant>.csv/.md
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

SUMMARY_KEYS = [
    ("points", "points"),
    ("clusters_raw", "candidate clusters (raw)"),
    ("cluster_points_raw", "clustered points (raw)"),
    ("clusters_adjusted", "candidate clusters (elevation-adjusted)"),
    ("cluster_points_adjusted", "clustered points (elevation-adjusted)"),
    ("share_of_raw_points_retained", "share of raw cluster points retained"),
    ("cluster_points_within_250m_raw", "headworks cluster points (raw)"),
    ("cluster_points_within_250m_adjusted", "headworks cluster points (adjusted)"),
    ("gradient_reported_velocity_mm_yr_per_km", "velocity-elevation gradient (mm/yr per km)"),
    ("spearman_reported_velocity", "Spearman rho velocity-elevation"),
    ("kappa_trend_mm_per_km_per_yr", "kappa trend (mm/km per yr)"),
    ("headworks_velocity_anomaly_raw_mm_yr", "headworks anomaly, raw (mm/yr)"),
    ("headworks_velocity_anomaly_adjusted_mm_yr", "headworks anomaly, adjusted (mm/yr)"),
]


def worker(settings_file: str, agent_dir: str, out_file: str) -> None:
    os.environ["KALI_SETTINGS"] = settings_file
    sys.path.insert(0, agent_dir)
    import numpy as np
    from kali_agent.config import settings
    from kali_agent.elevation import summary
    from kali_agent.screening import build
    from kali_agent.sites import zone_table

    s = settings()
    res = {"settings": settings_file, "products": {}}
    for key, p in s.products.items():
        if not p.path.is_file():
            res["products"][key] = {"missing": str(p.path)}
            continue
        for corr in ("none", "elevation"):
            build(key, corr, force=False)
        zones = {}
        for corr in ("none", "elevation"):
            try:
                tab, _ = zone_table(key, corr)
                zones[corr] = tab.to_dict(orient="records")
            except ValueError as exc:
                zones[corr] = [{"error": str(exc)}]
        res["products"][key] = {"summary": summary(key), "zones": zones}
        print(f"    {key}: done", flush=True)

    def clean(o):
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, list):
            return [clean(v) for v in o]
        if isinstance(o, (np.floating, float)):
            return None if not np.isfinite(o) else float(o)
        if isinstance(o, np.integer):
            return int(o)
        return o

    Path(out_file).write_text(json.dumps(clean(res), indent=1))


def _fmt(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        return f"{v:.3f}" if abs(v) < 1 else f"{v:.2f}"
    return str(v)


def side_by_side(results: dict[str, dict], variant: str, out_dir: Path) -> None:
    runs = list(results)
    rows = []
    products = []
    for r in results.values():
        for k in r.get("products", {}):
            if k not in products:
                products.append(k)
    for key in products:
        for skey, label in SUMMARY_KEYS:
            row = {"product": key, "quantity": label}
            for run in runs:
                p = results[run].get("products", {}).get(key, {})
                row[run] = p.get("summary", {}).get(skey)
            rows.append(row)
        zone_names = []
        for run in runs:
            for z in results[run].get("products", {}).get(key, {}).get("zones", {}).get("none", []):
                if z.get("zone") and z["zone"] not in zone_names:
                    zone_names.append(z["zone"])
        for corr, tag in (("none", "raw"), ("elevation", "adjusted")):
            for zn in zone_names:
                row = {"product": key, "quantity": f"{zn} rate, {tag} (mm/yr ± 95%)"}
                for run in runs:
                    zs = results[run].get("products", {}).get(key, {}).get("zones", {}).get(corr, [])
                    z = next((z for z in zs if z.get("zone") == zn), None)
                    if z is None:
                        row[run] = None
                    elif z.get("rate_mm_yr") is not None:
                        ci = z.get("ci95_mm_yr")
                        row[run] = (f"{z['rate_mm_yr']:+.2f} ± {ci:.2f}" if ci is not None else
                                    f"{z['rate_mm_yr']:+.2f}") + f" (n={z.get('points')})"
                    else:
                        row[run] = f"{z.get('note', '')} (n={z.get('points')})"
                rows.append(row)
    import csv
    csv_path = out_dir / f"answers_side_by_side_{variant}.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["product", "quantity"] + runs)
        w.writeheader()
        w.writerows(rows)
    md = [f"# Agent answers, SARscape vs open source (variant '{variant}')", "",
          "Columns: **sarscape** = products as used in the manuscript; **sarscape_matched** = the SARscape "
          "points that also have an open-source value, on the common dates; **opensource_at_sarscape** = "
          "open-source histories at exactly those points (same positions, elevations and zones); "
          "**opensource_native** = the open-source chain's own points. The like-for-like comparison is "
          "sarscape_matched vs opensource_at_sarscape.", "",
          "| product | quantity | " + " | ".join(runs) + " |", "|---|---|" + "---|" * len(runs)]
    for r in rows:
        md.append(f"| {r['product']} | {r['quantity']} | " + " | ".join(_fmt(r[run]) for run in runs) + " |")
    (out_dir / f"answers_side_by_side_{variant}.md").write_text("\n".join(md) + "\n", encoding="utf-8")
    print(f"  {csv_path}")


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--agent-dir", required=True, help="folder that contains kali_agent and settings.toml")
    ap.add_argument("--settings-dir", required=True, help="<export_dir>/agent_settings")
    ap.add_argument("--variant", default="no", help="tropospheric variant of the open-source run")
    ap.add_argument("--out", default=None, help="output folder (default: <settings-dir>/../answers)")
    ap.add_argument("--worker", nargs=2, metavar=("SETTINGS", "OUT"), help=argparse.SUPPRESS)
    a = ap.parse_args(argv)
    if a.worker:
        worker(a.worker[0], a.agent_dir, a.worker[1])
        return 0
    sdir = Path(a.settings_dir)
    out = Path(a.out) if a.out else sdir.parent / "answers"
    out.mkdir(parents=True, exist_ok=True)
    names = ["sarscape", f"sarscape_matched_{a.variant}", f"opensource_at_sarscape_{a.variant}",
             f"opensource_native_{a.variant}"]
    results = {}
    for name in names:
        f = sdir / f"{name}.toml"
        if not f.is_file():
            print(f"  skip {name}: {f} not found")
            continue
        target = out / f"answers_{name}.json"
        print(f"  {name} ...", flush=True)
        rc = subprocess.run([sys.executable, __file__, "--agent-dir", a.agent_dir, "--settings-dir", a.settings_dir,
                             "--worker", str(f), str(target)]).returncode
        if rc != 0 or not target.is_file():
            print(f"  {name}: FAILED (exit {rc})")
            continue
        results[name.replace(f"_{a.variant}", "")] = json.loads(target.read_text())
    if not results:
        print("No run succeeded.")
        return 1
    side_by_side(results, a.variant, out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
