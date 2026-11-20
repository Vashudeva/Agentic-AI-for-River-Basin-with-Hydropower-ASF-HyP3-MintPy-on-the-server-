# kgproc: open-source re-processing of the Kali Gandaki A analysis window

Two routes, one set of outputs (points in the agent's format, comparison with SARscape E-SBAS,
the agent's answers side by side):

| route | where | how | guide |
|---|---|---|---|
| **hyp3** | Windows (GEOK_Lab1), no admin rights | ASF HyP3 burst interferograms (ISCE2 run by ASF, ~20 m) + MintPy | [README_HYP3_WINDOWS.md](README_HYP3_WINDOWS.md) |
| **isce** | Linux / WSL2 | ISCE2 + MiaplPy (PS + DS at full resolution) + MintPy, from your own SAFE files | this file |

This folder re-processes the **same Sentinel-1 images** that went into the SARscape E-SBAS
run, without ENVI/SARscape, and compares the two results point by point and through the
agent's own answers.

```
Sentinel-1 SAFE (.zip)  ──► ISCE2 topsStack ──► MiaplPy ──────────────► MintPy ─────────► kgproc products ──► compare
  same 64 + 63 dates        coregistered        PS + DS phase linking,   corrections,       agent-format         statistics,
  (read from CG_report)     full-res SLCs,      single-reference         geocoding to       GeoPackages          figures, report
                            bursts of the       network on the SARscape  SARscape's grid    + agent settings  ──► answers.py (agent,
                            window only         reference date, SNAPHU   (0.000125°)                              Windows): Table 5,
                                                                                                                   zones A–F, side by side
```

Everything is Python (ISCE2's heavy parts are compiled C/C++/Fortran called from Python).
No language model takes part in the processing; the agent only reads the finished products.

## What is matched to the SARscape run, and what is not

| | SARscape E-SBAS (manuscript) | kgproc |
|---|---|---|
| Images | 64 asc (track 158, IW2) / 63 dsc (track 92, IW1) | the same dates, read from `CG_report.txt` and checked |
| Window | crop of each stack | same `window_snwe` per track |
| Network | one reference image (10 Oct 2022 / 12 Sep 2022) | MiaplPy `single_reference` on the same dates |
| PS / DS | E-SBAS (SARscape's own DS handling) | MiaplPy phase linking (sequential EMI, KS test for neighbours) |
| Unwrapping | SARscape | SNAPHU |
| Atmosphere | SARscape atmospheric filtering (`slant_atm_meta` in the inversion folder) | variant `no` (none), `pyaps` (ERA5), `gacos`, `height_correlation` |
| Reference point | separate: asc 1,593 m, dsc 739 m | **one** for both tracks (needed for the decomposition) |
| DEM | [DEM placeholder in the manuscript] | Copernicus GLO-30, or your DEM via `isce.dem_file` |
| Geocoding grid | 0.000125° | 0.000125°, same extent for both tracks |
| Decomposition | d_N = 0, each LOS series interpolated to the other's dates (124 dates) | the same (checked: 124 dates) |

So the comparison measures the effect of the **processor** (coregistration, DS estimation,
unwrapping, atmosphere handling, reference), not of the data.

## 1. One-time setup on GEOK_Lab1 (Windows → WSL2)

ISCE2 and MiaplPy run only on Linux. On Windows this means WSL2. **Turning WSL on is a
system change (admin rights, reboot) and lies outside the G:\InSAR\Kaligandaki\Langgraph rule**;
everything after that can live inside the folder.

In an administrator PowerShell:

```powershell
wsl --install --no-distribution            # enables WSL2; reboot afterwards
wsl --install -d Ubuntu-24.04              # temporary copy on C:
wsl --export Ubuntu-24.04 G:\InSAR\Kaligandaki\Langgraph\wsl\ubuntu.tar
wsl --unregister Ubuntu-24.04
wsl --import KaliInSAR G:\InSAR\Kaligandaki\Langgraph\wsl\KaliInSAR G:\InSAR\Kaligandaki\Langgraph\wsl\ubuntu.tar
del G:\InSAR\Kaligandaki\Langgraph\wsl\ubuntu.tar
```

The Linux file system then is `G:\InSAR\Kaligandaki\Langgraph\wsl\KaliInSAR\ext4.vhdx`
(it grows to roughly 300 GB during processing). Optional: give WSL more memory in
`C:\Users\<you>\.wslconfig` (`[wsl2]` / `memory=48GB` / `processors=12`); 32 GB or more is advisable.

Copy this folder to `G:\InSAR\Kaligandaki\Langgraph\opensource_insar`, then in WSL
(`wsl -d KaliInSAR`):

```bash
sudo apt-get update && sudo apt-get install -y curl git bzip2      # if missing
cd /mnt/g/InSAR/Kaligandaki/Langgraph/opensource_insar
bash setup_env.sh                       # conda environment ~/kali-insar-env (10-20 min)
```

For the ERA5 variant create `~/.cdsapirc` with your Copernicus Climate Data Store key
(format on the CDS "how to API" page, currently `url: https://cds.climate.copernicus.eu/api`
and `key: <personal access token>`), and accept the ERA5 licence on the CDS website once.

## 2. Configure

Edit `config.toml`. Normally only:

* `tracks.*.slc_dir` – the folder with the original `S1*_IW_SLC__1S*.zip` files
  (both tracks may share one folder; sub-folders are searched). The SARscape-imported
  `*_slc_list` files cannot be used: ISCE2 needs the original SAFE products.
* `paths.work_dir` – keep it on the Linux side (default `~/kali_work`): ISCE2 on `/mnt/g` is very slow.
* `mintpy.reference_lalo` – a coherent, stable pixel in **both** windows (default: SARscape's descending reference).
* `mintpy.troposphere` – the variants to produce (`["no", "pyaps"]` by default).
* `isce.dem_file` – set it to the DEM SARscape used if you want to remove the DEM as a difference.

## 3. Run

```bash
cd /mnt/g/InSAR/Kaligandaki/Langgraph/opensource_insar
PY=~/kali-insar-env/bin/python

$PY -m kgproc check                  # environment, inputs, free space
$PY -m kgproc inventory              # finds the 64 + 63 scenes, stops if any is missing
$PY -m kgproc aux                    # precise orbits (ASF, no account) + DEM (checked: a failed download is refused)
$PY -m kgproc isce                   # ISCE2 stack, resumable run file by run file
$PY -m kgproc miaplpy                # phase linking -> unwrapping -> network inversion
$PY -m kgproc correct                # MintPy corrections + geocoding, once per variant
$PY -m kgproc products               # point products + agent settings
$PY -m kgproc compare                # statistics, figures, report
```

or all at once: `nohup $PY -m kgproc all > kgproc.log 2>&1 &`. Every step can be rerun;
finished parts are skipped (markers in `<work_dir>/_done`). `--track ascending` or
`--variant no` restrict a step. Logs: `<work_dir>/<track>/logs/`.

Rough cost per track on an 8–12-core workstation: ISCE2 several hours; MiaplPy phase linking
at full resolution about one day (it scales with window size × number of images);
about 150 GB of disk. These are estimates, not measurements.

## 4. The agent's answers on both product sets (Windows)

With the agent's own Python (no WSL needed), in `cmd`:

```bat
cd /d G:\InSAR\Kaligandaki\Langgraph\opensource_insar
G:\InSAR\Kaligandaki\Langgraph\.venv\Scripts\python.exe kgproc\answers.py ^
   --agent-dir G:\InSAR\Kaligandaki\Langgraph\insar-hydro-agent\agent ^
   --settings-dir G:\InSAR\Kaligandaki\Langgraph\opensource\agent_settings ^
   --variant no
```

It runs the agent's `build` (Table 5 numbers, elevation check, zones A–F double
differences) four times, each with its own settings file and work folder:

| run | products |
|---|---|
| `sarscape` | SARscape files as in the manuscript (should reproduce Table 5) |
| `sarscape_matched` | SARscape points that also have an open-source value, common dates |
| `opensource_at_sarscape` | open-source histories at exactly those points (same positions, elevations, zones) |
| `opensource_native` | the open-source chain's own points |

`sarscape_matched` vs `opensource_at_sarscape` is the like-for-like comparison.
Result: `opensource\answers\answers_side_by_side_<variant>.md/.csv`. Run `compare` again
afterwards so the report links to it. The live agent, its cache and the web page are not touched.

## 5. Outputs (`G:\InSAR\Kaligandaki\Langgraph\opensource`)

```
points\<variant>\os_native_<product>.gpkg          open-source points (QGIS / agent)
points\<variant>\os_at_sarscape_<product>.gpkg     open-source values at the SARscape points
points\<variant>\sarscape_matched_<product>.gpkg   the matching SARscape subset
points\<variant>\offset_scan.csv                   geolocation check
comparison\<variant>\comparison_report.md          tables + figures
comparison\<variant>\velocity_agreement.csv ...    all tables as CSV
agent_settings\*.toml                              settings for the answers runs
answers\answers_side_by_side_<variant>.md          the agent's numbers, side by side
```

All products: displacements in mm re-zeroed at the first date; LOS positive towards the
satellite, vertical positive up, east positive east; `coherence` = MintPy temporal coherence;
`scatterer` = PS/DS from MiaplPy's PS mask.

## 6. Reading the comparison

* **Velocity agreement** (section 2 of the report): median difference = offset between
  the chains after common re-referencing; NMAD = point-level scatter. Set it against the
  theoretical precision in Table 4 (0.56 / 0.55 mm/yr).
* **Per-date difference** (lower-right panels): a pattern common to all points and dates
  points to atmosphere, reference handling or unwrapping, not to ground motion.
* **Elevation gradient** (section 4): the direct test of the manuscript's Section 4.2. If
  the height-dependent bias (dsc ≈ −10, vertical ≈ −12 mm/yr per km) is residual
  stratified troposphere, it should shrink in the `pyaps`/`gacos` variants.
* **Offset scan** (section 5): which positions of the decomposed products (lon/lat fields
  or the shifted geometry of the manuscript copies) agree with independent geocoding.
* **Agent answers**: whether zone D (−3.3 mm/yr) and the cluster counts survive a change
  of processor.

## 7. Self-check (minutes, no Sentinel-1 data needed)

```bash
$PY tests/synthetic_from_sarscape.py          # writes a variant 'synthetic' from the SARscape LOS points
$PY -m kgproc products --variant synthetic
$PY -m kgproc compare  --variant synthetic
```

Tested on 20 Nov 2026 with the SARscape products: LOS differences zero except where two
SARscape points fall into one grid cell; kgproc's decomposition of the SARscape LOS data
against SARscape's own vertical/east: r = 0.96 / 0.94, NMAD 0.33 / 0.42 mm/yr, slope 0.99;
the agent run `sarscape` reproduced Table 5 (93/1,601, 97/1,605, 113/2,102, 97/1,662;
vertical 38/779 → 35/527) and zone D vertical −3.32 ± 0.45 mm/yr.
Not tested here: the ISCE2, MiaplPy and MintPy runs themselves (no SAFE files in the test);
their command lines and option names were checked against isce2 2.6.5, MiaplPy 0.2.1 and
MintPy 1.6.4.

## 8. If this goes into the manuscript

Software to cite (verify each record, e.g. via Crossref, before inserting with Mendeley Cite):
ISCE2 topsStack (Fattahi, Agram and Simons, 2017, IEEE TGRS, network-based enhanced spectral
diversity); MiaplPy (Mirzaee, Amelung and Fattahi, 2023, Computers & Geosciences 171, 105291);
MintPy (Yunjun, Fattahi and Amelung, 2019, Computers & Geosciences 133, 104331); PyAPS/ERA5
(Jolivet et al., 2011, already in the reference list); SNAPHU (Chen and Zebker).
