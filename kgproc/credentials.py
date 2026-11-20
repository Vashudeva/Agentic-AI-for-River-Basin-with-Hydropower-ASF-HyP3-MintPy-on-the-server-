"""Credentials kept inside work_dir (never printed).

NASA Earthdata Login token (for ASF HyP3): `python -m kgproc set-earthdata-token`
  stored in <work_dir>/env/earthdata_token.dpapi, encrypted with Windows DPAPI for the
  current Windows user (plain file with owner-only rights on Linux). The environment
  variable EARTHDATA_TOKEN, if set, is used instead.
Copernicus CDS key (ERA5 for the 'pyaps' variant): `python -m kgproc set-cds-key`
  written to PyAPS's own model.cfg inside the Python environment, so nothing goes to the
  user profile. An existing ~/.cdsapirc takes precedence in PyAPS.
"""
from __future__ import annotations

import configparser
import getpass
import os
from pathlib import Path

from .config import Config


def _dpapi(data: bytes, protect: bool) -> bytes:
    import ctypes
    from ctypes import wintypes

    class BLOB(ctypes.Structure):
        _fields_ = [("cbData", wintypes.DWORD), ("pbData", ctypes.POINTER(ctypes.c_char))]

    buf = ctypes.create_string_buffer(data, len(data))
    blob_in = BLOB(len(data), ctypes.cast(buf, ctypes.POINTER(ctypes.c_char)))
    blob_out = BLOB()
    fn = ctypes.windll.crypt32.CryptProtectData if protect else ctypes.windll.crypt32.CryptUnprotectData
    if not fn(ctypes.byref(blob_in), None, None, None, None, 0, ctypes.byref(blob_out)):
        raise OSError("Windows DPAPI call failed")
    try:
        return ctypes.string_at(blob_out.pbData, blob_out.cbData)
    finally:
        ctypes.windll.kernel32.LocalFree(ctypes.cast(blob_out.pbData, ctypes.c_void_p))


def _token_file(cfg: Config) -> Path:
    return cfg.work_dir / "env" / ("earthdata_token.dpapi" if os.name == "nt" else "earthdata_token")


def set_token(cfg: Config) -> None:
    print("Create a token at https://urs.earthdata.nasa.gov (Generate Token); it is valid for 60 days.")
    tok = getpass.getpass("NASA Earthdata Login token (input hidden): ").strip()
    if not tok:
        raise SystemExit("No token entered.")
    f = _token_file(cfg)
    f.parent.mkdir(parents=True, exist_ok=True)
    if os.name == "nt":
        f.write_bytes(_dpapi(tok.encode(), protect=True))
    else:
        f.write_text(tok)
        os.chmod(f, 0o600)
    print(f"Token stored ({'DPAPI-encrypted, ' if os.name == 'nt' else ''}{f}).")


def load_token(cfg: Config) -> str | None:
    if os.environ.get("EARTHDATA_TOKEN"):
        return os.environ["EARTHDATA_TOKEN"].strip()
    f = _token_file(cfg)
    if not f.is_file():
        return None
    return (_dpapi(f.read_bytes(), protect=False).decode() if os.name == "nt" else f.read_text()).strip()


def set_cds_key() -> None:
    import pyaps3
    cfg_file = Path(pyaps3.__file__).with_name("model.cfg")
    print("Copernicus Climate Data Store: copy your Personal Access Token from your CDS profile page,")
    print("and accept the ERA5 (pressure levels) licence on the CDS website once.")
    key = getpass.getpass("CDS personal access token (input hidden): ").strip()
    if not key:
        raise SystemExit("No key entered.")
    cp = configparser.RawConfigParser(delimiters="=")
    cp.read(cfg_file)
    if not cp.has_section("CDS"):
        cp.add_section("CDS")
    cp.set("CDS", "key", key)
    with open(cfg_file, "w") as fh:
        cp.write(fh)
    print(f"CDS key written to {cfg_file}")
