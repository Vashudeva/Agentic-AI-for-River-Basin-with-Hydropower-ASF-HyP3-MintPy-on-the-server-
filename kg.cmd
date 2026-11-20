@echo off
rem kgproc launcher, HyP3 route on Windows. Usage:  kg <step> [options]   e.g.  kg check
rem Everything stays inside G:\InSAR\Kaligandaki\Langgraph.
setlocal
set "ROOT=G:\InSAR\Kaligandaki\Langgraph"
set "MAMBA_ROOT_PREFIX=%ROOT%\micromamba"
set "KGENV=%ROOT%\envs\kali-hyp3"
set "KGWORK=%ROOT%\opensource_work"
set "TEMP=%KGWORK%\tmp"
set "TMP=%KGWORK%\tmp"
set "MPLCONFIGDIR=%KGWORK%\env\matplotlib"
set "PROJ_USER_WRITABLE_DIRECTORY=%KGWORK%\env\proj"
set "PYTHONUNBUFFERED=1"
if not exist "%TEMP%" mkdir "%TEMP%"
if not exist "%MAMBA_ROOT_PREFIX%\micromamba.exe" (
  echo micromamba not found: run setup_windows.cmd first.
  exit /b 1
)
cd /d "%~dp0"
"%MAMBA_ROOT_PREFIX%\micromamba.exe" run -p "%KGENV%" python -m kgproc -c config_hyp3.toml %*
exit /b %ERRORLEVEL%
