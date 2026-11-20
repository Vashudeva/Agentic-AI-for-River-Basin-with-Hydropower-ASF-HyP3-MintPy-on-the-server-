"""kgproc: open-source re-processing of the Kali Gandaki A analysis window and its
comparison with the SARscape E-SBAS products.

Steps (python -m kgproc <step>):
  check      environment and inputs
  inventory  find the Sentinel-1 files of the SARscape dates, link them per track
  aux        precise orbits and DEM
  isce       ISCE2 topsStack: coregistered full-resolution SLC stack of the window
  miaplpy    MiaplPy: phase linking, interferograms, unwrapping, network inversion
  correct    MintPy corrections and geocoding, one run per tropospheric variant
  products   decomposition and point products in the agent's format
  compare    SARscape vs open-source statistics, figures and report
  all        every step in order

The answers of the agent on both product sets come from kgproc/answers.py, which runs
with the agent's own Python on Windows (see README.md).

Light modules (config, compare, answers) import only numpy/pandas/geopandas, so that
the comparison also runs in the agent's Windows environment.
"""
__version__ = "1.0.0"
