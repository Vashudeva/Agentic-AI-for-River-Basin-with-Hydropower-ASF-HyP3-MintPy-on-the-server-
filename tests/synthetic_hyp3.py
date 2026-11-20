"""Dry run of the HyP3 route without HyP3 credits or downloads.

Builds HyP3-like INSAR_ISCE_MULTI_BURST product zips (UTM GeoTIFFs + metadata .txt) whose
unwrapped phase is computed from the SARscape LOS histories, then runs the real code of
the route on them: unpack + clip, MintPy (load, inversion, corrections), resampling to the
SARscape grid, points, comparison. Everything goes to separate 'selftest' folders.

    kg check
    %ROOT%\\micromamba\\micromamba.exe run -p %ROOT%\\envs\\kali-hyp3 python tests\\synthetic_hyp3.py -c config_hyp3.toml

Expected: line-of-sight velocities equal to SARscape up to the 40 m test pixels (the test
uses 40 m to stay small); vertical/east close to SARscape's own decomposition.
This checks that MintPy and the whole route run on this machine before credits are spent.
"""
from __future__ import annotations

import argparse
import shutil
import sys
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from kgproc import compare, hyp3_run, products               # noqa: E402
from kgproc.config import TRACKS, load                       # noqa: E402
from kgproc.inventory import parse_cg_report                 # noqa: E402

WAVELENGTH = 0.05546576
HEADING = {"ascending": 347.8, "descending": 192.2}


def write_tif(path: Path, arr: np.ndarray, gt, epsg: int, dtype) -> None:
    from osgeo import gdal, osr
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(str(path), arr.shape[1], arr.shape[0], 1, dtype)
    ds.SetGeoTransform(gt)
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(epsg)
    ds.SetProjection(srs.ExportToWkt())
    ds.GetRasterBand(1).WriteArray(arr)
    ds = None


def main(argv=None) -> int:
    from pyproj import Transformer
    from osgeo import gdal
    from scipy.spatial import cKDTree
    ap = argparse.ArgumentParser()
    ap.add_argument("-c", "--config", default=None)
    ap.add_argument("--res", type=float, default=40.0, help="test pixel size in m")
    ap.add_argument("--max-dates", type=int, default=0, help="use only the first N dates (quicker)")
    a = ap.parse_args(argv)
    cfg = load(a.config)
    cfg.work_dir = cfg.work_dir / "selftest"
    cfg.export_dir = cfg.export_dir / "selftest"
    cfg.hyp3 = dict(cfg.hyp3, looks="10x2")
    rng = np.random.default_rng(1)
    epsg = 32644

    for key in TRACKS:
        t = cfg.tracks[key]
        fr, dcols, lon, lat, _ = products.read_sarscape(cfg.sarscape_products[key])
        if a.max_dates:
            dcols = dcols[:a.max_dates]
        dates = [c[2:] for c in dcols]
        D = fr[dcols].apply(pd.to_numeric, errors="coerce").to_numpy(float) / 1000.0
        S, N, W, E = t.window_snwe
        tr = Transformer.from_crs("EPSG:4326", f"EPSG:{epsg}", always_xy=True)
        xs, ys = tr.transform(np.array([W, E, W, E]), np.array([S, S, N, N]))
        x0, y1 = np.floor(min(xs) / a.res) * a.res - 4 * a.res, np.ceil(max(ys) / a.res) * a.res + 4 * a.res
        cols = int((np.ceil(max(xs) / a.res) * a.res + 4 * a.res - x0) / a.res)
        rows = int((y1 - (np.floor(min(ys) / a.res) * a.res - 4 * a.res)) / a.res)
        gt = (x0, a.res, 0.0, y1, 0.0, -a.res)
        px, py = tr.transform(lon, lat)
        X, Y = np.meshgrid(x0 + (np.arange(cols) + 0.5) * a.res, y1 - (np.arange(rows) + 0.5) * a.res)
        dist, idx = cKDTree(np.c_[px, py]).query(np.c_[X.ravel(), Y.ravel()], distance_upper_bound=0.75 * a.res)
        valid = np.isfinite(dist).reshape(rows, cols)
        idx = np.where(np.isfinite(dist), idx, 0).reshape(rows, cols)
        series = D[idx]                                                     # rows, cols, dates
        num = lambda f: pd.to_numeric(fr[f], errors="coerce").to_numpy(float)
        h = np.where(valid, num("Z")[idx], np.nanmedian(num("Z"))).astype(np.float32)
        inc = np.where(valid, num("ILOS")[idx], np.nanmedian(num("ILOS")))
        alos = np.where(valid, num("ALOS")[idx], np.nanmedian(num("ALOS")))
        az = np.sign(alos) * (180.0 - np.abs(alos))                       # MintPy LOS azimuth
        lv_theta = np.radians(90.0 - inc).astype(np.float32)               # HyP3: from horizontal, radians
        lv_phi = np.radians(az + 90.0).astype(np.float32)                  # HyP3: from east, radians
        pairs = hyp3_run.network(dates)
        bid = "000000s1n00-336254s2n02-000000s3n00" if t.swath == 2 else "141234s1n02-000000s2n00-000000s3n00"
        zdir = cfg.hyp3_dir(key) / "_synthetic_zips"
        zdir.mkdir(parents=True, exist_ok=True)
        print(f"{key}: {rows} x {cols} pixels of {a.res:.0f} m, {int(valid.sum())} with SARscape data, "
              f"{len(dates)} dates, {len(pairs)} pairs")
        k_of = {d: k for k, d in enumerate(dates)}
        for n, (d1, d2) in enumerate(pairs):
            product = f"S1_{t.relative_orbit:03d}_{bid}_IW_{d1}_{d2}_VV_INT40_{n:04X}"
            work = zdir / product
            work.mkdir(exist_ok=True)
            dphase = -(4 * np.pi / WAVELENGTH) * (series[..., k_of[d2]] - series[..., k_of[d1]])
            unw = np.where(valid, dphase, rng.uniform(-np.pi, np.pi, (rows, cols)) * 3).astype(np.float32)
            corr = np.where(valid, 0.9, 0.15).astype(np.float32)
            write_tif(work / f"{product}_unw_phase.tif", unw, gt, epsg, gdal.GDT_Float32)
            write_tif(work / f"{product}_corr.tif", corr, gt, epsg, gdal.GDT_Float32)
            write_tif(work / f"{product}_conncomp.tif", np.ones((rows, cols), np.uint8), gt, epsg, gdal.GDT_Byte)
            if n == 0:
                write_tif(work / f"{product}_dem.tif", h, gt, epsg, gdal.GDT_Float32)
                write_tif(work / f"{product}_lv_theta.tif", lv_theta, gt, epsg, gdal.GDT_Float32)
                write_tif(work / f"{product}_lv_phi.tif", lv_phi, gt, epsg, gdal.GDT_Float32)
                write_tif(work / f"{product}_water_mask.tif", np.ones((rows, cols), np.uint8), gt, epsg, gdal.GDT_Byte)
            (work / f"{product}.txt").write_text(
                f"Reference Granule: synthetic\nSecondary Granule: synthetic\nBaseline: {rng.uniform(-150, 150):.3f}\n"
                f"UTC time: 45000.0\nHeading: {HEADING[key]}\nSpacecraft height: 693000.0\n"
                f"Earth radius at nadir: 6337286.6\nSlant range near: 800000.0\nRange looks: 10\n"
                f"Azimuth looks: 2\nUnwrapping type: snaphu_mcf\n")
            zp = zdir / f"{product}.zip"
            with zipfile.ZipFile(zp, "w") as z:
                for f in work.iterdir():
                    z.write(f, f"{product}/{f.name}")
            shutil.rmtree(work)
            hyp3_run._unpack(cfg, key, zp, need_geometry=(n == 0))
            zp.unlink()
        shutil.rmtree(zdir, ignore_errors=True)
        hyp3_run.run_mintpy(cfg, key, "no")

    products.run_step(cfg, "hyp3_no")
    compare.run_step(cfg, "hyp3_no")
    print(f"\nSelf-test report: {cfg.report_dir('hyp3_no') / 'comparison_report.md'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
