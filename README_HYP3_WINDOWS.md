# HyP3 route on Windows (GEOK_Lab1, no admin rights, no Linux)

```
SARscape dates (CG_report.txt)        ASF HyP3 (ISCE2 on ASF's computers)       GEOK_Lab1 (Windows)
+ SARscape windows             ──►    INSAR_ISCE_MULTI_BURST, 5x1 looks ~20 m ──► MintPy ──► points ──► comparison ──► agent answers
kg hyp3-search / hyp3-submit          small-baseline network                     kg hyp3-download / hyp3-mintpy / products / compare
```

* Same acquisitions as the ENVI run: the dates are read from SARscape's connection-graph
  reports; the bursts are those of IW2 (ascending, track 158) and IW1 (descending, track 92)
  that cover the SARscape crop windows. HyP3 uses ASF's copies of these images, not the
  files on F:.
* Network: each date with its 3 next dates plus one pair of about one year per date:
  239 pairs ascending + 236 descending = 475 jobs (186 + 183 without the yearly pairs,
  `annual_pairs = false`).
* Everything is installed and written inside `G:\InSAR\Kaligandaki\Langgraph`
  (environment `envs\kali-hyp3`, work `opensource_work`, results `opensource`).

## Before you start

* A free NASA Earthdata Login account (urs.earthdata.nasa.gov) and a token from its
  "Generate Token" page (valid 60 days).
* Internet access from the server to `*.asf.alaska.edu`, `cmr.earthdata.nasa.gov`,
  `urs.earthdata.nasa.gov`, `conda.anaconda.org`, `micro.mamba.pm` (and the CDS for ERA5).
* HyP3 has a monthly credit allowance. `kg hyp3-submit` prints the estimated cost and your
  remaining credits and submits nothing unless you add `--yes`.
* ASF policy: **HyP3 jobs are publicly visible**, and finished products are deleted after
  about two weeks, so download them promptly (expired ones can be resubmitted).

## 1. Setup (once)

```bat
cd /d G:\InSAR\Kaligandaki\Langgraph
git clone https://github.com/Vashudeva/Agentic-AI-for-River-Basin-with-Hydropower-ASF-HyP3-MintPy-on-the-server-.git hyp3-mintpy
cd hyp3-mintpy
setup_windows.cmd                 &rem micromamba + Python environment (~500 MB, 10-20 min)
kg set-earthdata-token            &rem stored DPAPI-encrypted in opensource_work\env
kg set-cds-key                    &rem only for the ERA5 variant (written into the environment)
kg check
```

(Set `$env:GIT_CONFIG_GLOBAL` / `set GIT_CONFIG_GLOBAL=G:\InSAR\Kaligandaki\Langgraph\work\env\gitconfig`
first, as for the agent.) Paths are in `config_hyp3.toml` and at the top of `kg.cmd` and
`setup_windows.cmd`; change them there if the folder ever moves.

## 2. Dry run without credits (recommended, ~15 min)

```bat
G:\InSAR\Kaligandaki\Langgraph\micromamba\micromamba.exe run -p G:\InSAR\Kaligandaki\Langgraph\envs\kali-hyp3 python tests\synthetic_hyp3.py -c config_hyp3.toml
```

It writes HyP3-format interferograms computed from the SARscape histories (40 m test pixels)
and runs the real code on them: unpack and clip, MintPy, resampling, points, comparison.
Output in `opensource_work\selftest` and `opensource\selftest`. If this finishes, MintPy
works on the server. (Tested on 20 Nov 2026 in Linux with the same package versions:
line-of-sight velocities r = 0.91-0.92, vertical 0.93 against SARscape, the residual
coming from the coarse 40 m test pixels.)

## 3. Run

```bat
kg hyp3-search                    &rem bursts per date; stops if dates are missing
kg hyp3-submit                    &rem shows pairs, cost and credits; nothing submitted
kg hyp3-submit --yes              &rem submits (rerunning it never duplicates jobs; failed ones are resubmitted)
kg hyp3-status                    &rem progress
kg hyp3-download                  &rem downloads finished products, keeps only the window
kg hyp3-download --wait           &rem ... or waits for the rest (up to wait_hours)
kg hyp3-mintpy                    &rem MintPy, variants "no" and "pyaps" (ERA5)
kg products                       &rem points: native, at the SARscape points, agent settings
kg compare                        &rem report: opensource\comparison\hyp3_<variant>\comparison_report.md
```

`--track ascending` or `--variant hyp3_no` restrict a step. Logs: `opensource_work\<track>\logs`.

The agent's answers (Table 5 numbers, zones A–F) on both product sets, with the agent's Python:

```bat
G:\InSAR\Kaligandaki\Langgraph\.venv\Scripts\python.exe kgproc\answers.py ^
  --agent-dir G:\InSAR\Kaligandaki\Langgraph\insar-hydro-agent\agent ^
  --settings-dir G:\InSAR\Kaligandaki\Langgraph\opensource\agent_settings --variant hyp3_no
```

Then `kg compare` again so the report links the answers. The live agent, its cache and the
web page are not touched.

## 4. What to expect

* Disk: a few GB (clipped interferograms ~3 GB, MintPy stacks ~3 GB per track); each product
  zip is deleted after the window is cut out (`keep_zips = false`).
* Time: HyP3 processing depends on ASF's queue; MintPy on the server takes minutes to about
  an hour per track. These are estimates.
* Resolution matters for the headworks zones: zone D (0.7 ha) holds about 17 pixels at 20 m.
  In the dry run at 40 m it held only 4 and its rate was diluted (descending −3.3 → −1.8 mm/yr)
  by resampling alone, so read zone-level differences together with the pixel size.
* Differences from SARscape by design: small-baseline pairs instead of one reference image,
  multilooked pixels instead of E-SBAS PS/DS points, ISCE2 coregistration and SNAPHU
  unwrapping run by ASF, Copernicus DEM, one common reference point for both tracks,
  optional ERA5 correction.

## 5. Network problems

The server's agent already uses the Windows certificate store (`truststore`); kgproc does the
same for its own Python downloads. If the university inspects TLS traffic, Python tools can be
pointed to the university root certificate (exported to a .pem file inside the Langgraph folder):

```bat
set SSL_CERT_FILE=G:\InSAR\Kaligandaki\Langgraph\work\env\university_root.pem
set REQUESTS_CA_BUNDLE=%SSL_CERT_FILE%
```

If `setup_windows.cmd` itself (curl or micromamba) fails with certificate errors, keep the
error message: the fix depends on the proxy.

## 6. Citing (verify each record before inserting with Mendeley Cite)

HyP3 (ASF HyP3 software and the burst InSAR product guide), ISCE2 topsStack (Fattahi, Agram
and Simons, 2017), MintPy (Yunjun, Fattahi and Amelung, 2019, Computers & Geosciences 133,
104331), SNAPHU (Chen and Zebker), ERA5/PyAPS (Jolivet et al., 2011).
