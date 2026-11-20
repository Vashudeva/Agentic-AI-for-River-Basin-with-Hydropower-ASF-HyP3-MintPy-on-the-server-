"""Configuration (config.toml) and folder layout."""
from __future__ import annotations

import tomllib
from dataclasses import dataclass
from pathlib import Path, PurePosixPath, PureWindowsPath

TRACKS = ("ascending", "descending")
DECOMPOSED = ("vertical", "east")
PRODUCTS = TRACKS + DECOMPOSED


@dataclass(frozen=True)
class Track:
    key: str
    slc_dir: Path
    slc_recursive: bool
    relative_orbit: int
    pass_direction: str
    swath: int
    reference_date: str
    window_snwe: tuple[float, float, float, float]
    sarscape_dates: Path
    sarscape_product: Path


class Config:
    def __init__(self, path: str | Path):
        self.path = Path(path).resolve()
        with open(self.path, "rb") as fh:
            self.raw = raw = tomllib.load(fh)
        p = raw["paths"]
        self.work_dir = Path(p["work_dir"]).expanduser()
        self.export_dir = Path(p["export_dir"]).expanduser()
        self.sarscape_root = Path(p["sarscape_root"]).expanduser()
        self.agent_dir = Path(p["agent_dir"]).expanduser() if p.get("agent_dir") else None
        self.tracks: dict[str, Track] = {}
        for key in TRACKS:
            t = raw["tracks"][key]
            self.tracks[key] = Track(
                key=key,
                slc_dir=Path(t.get("slc_dir", "")).expanduser(),
                slc_recursive=bool(t.get("slc_recursive", True)),
                relative_orbit=int(t["relative_orbit"]),
                pass_direction=str(t["pass_direction"]).upper(),
                swath=int(t["swath"]),
                reference_date=str(t["reference_date"]),
                window_snwe=tuple(float(v) for v in t["window_snwe"]),
                sarscape_dates=self.sarscape(t["sarscape_dates"]),
                sarscape_product=self.sarscape(t["sarscape_product"]),
            )
        dec = raw.get("sarscape_decomposed", {})
        self.sarscape_products = {k: self.tracks[k].sarscape_product for k in TRACKS}
        for k in DECOMPOSED:
            if dec.get(k):
                self.sarscape_products[k] = self.sarscape(dec[k])
        self.isce = raw.get("isce", {})
        self.miaplpy = raw.get("miaplpy", {})
        self.mintpy = raw.get("mintpy", {})
        self.export = raw.get("export", {})
        self.compare = raw.get("compare", {})
        self.hyp3 = raw.get("hyp3", {})
        self.route = str(raw.get("route", "isce")).lower()      # "isce" (ISCE2 + MiaplPy) or "hyp3"

    # ------------------------------------------------------------------ paths
    def sarscape(self, rel: str) -> Path:
        """SARscape paths may be written with backslashes (as in the agent settings)."""
        p = Path(str(rel).replace("\\", "/"))
        return p if p.is_absolute() else self.sarscape_root / p

    def track_dir(self, key: str) -> Path:
        return self.work_dir / key

    def isce_dir(self, key: str) -> Path:
        return self.track_dir(key) / "isce"

    def miaplpy_dir(self, key: str) -> Path:
        return self.track_dir(key) / "miaplpy"

    def network_dir(self, key: str) -> Path:
        net = self.miaplpy.get("network", "single_reference")
        return self.miaplpy_dir(key) / f"network_{net}"

    def geo_dir(self, key: str, variant: str) -> Path:
        """Geocoded MintPy files of one tropospheric variant (copied after each run)."""
        return self.track_dir(key) / "geo" / variant

    def points_dir(self, variant: str) -> Path:
        return self.export_dir / "points" / variant

    def report_dir(self, variant: str) -> Path:
        return self.export_dir / "comparison" / variant

    def hyp3_dir(self, key: str) -> Path:
        """Downloaded and clipped HyP3 interferograms of one track."""
        return self.track_dir(key) / "hyp3"

    def hyp3_mintpy_dir(self, key: str) -> Path:
        return self.track_dir(key) / "mintpy_hyp3"

    @property
    def tropo(self) -> list[str]:
        v = self.mintpy.get("troposphere", ["no"])
        return [v] if isinstance(v, str) else list(v)

    @property
    def variants(self) -> list[str]:
        """Product variants: one per tropospheric correction; HyP3 ones are prefixed 'hyp3_'."""
        return [f"hyp3_{v}" for v in self.tropo] if self.route == "hyp3" else self.tropo

    def marker(self, key: str, step: str) -> Path:
        return self.work_dir / "_done" / f"{key}.{step}"


def load(path: str | Path | None = None) -> Config:
    if path is None:
        path = Path(__file__).resolve().parents[1] / "config.toml"
    return Config(path)


def local_path(path: str) -> Path:
    """A path written in a settings file, usable on this machine (Windows or WSL)."""
    import os
    import re
    s = str(path)
    if os.name != "nt" and re.match(r"^[A-Za-z]:[\\/]", s):
        return Path(to_posix(s))
    return Path(s)


def to_windows(path: Path | str) -> str:
    """/mnt/g/InSAR/x -> G:\\InSAR\\x (paths the Windows agent must open)."""
    p = PurePosixPath(str(path))
    parts = p.parts
    if len(parts) >= 3 and parts[0] == "/" and parts[1] == "mnt" and len(parts[2]) == 1:
        return str(PureWindowsPath(f"{parts[2].upper()}:\\", *parts[3:]))
    return str(path)


def to_posix(path: str) -> str:
    """G:\\InSAR\\x -> /mnt/g/InSAR/x (inverse of to_windows)."""
    s = str(path)
    if len(s) >= 2 and s[1] == ":":
        rest = s[2:].replace("\\", "/").lstrip("/")
        return f"/mnt/{s[0].lower()}/{rest}"
    return s.replace("\\", "/")
