@echo off
rem One-time setup of the HyP3 route on GEOK_Lab1 (Windows, no admin rights needed).
rem Installs micromamba and a Python environment inside G:\InSAR\Kaligandaki\Langgraph.
setlocal
set "ROOT=G:\InSAR\Kaligandaki\Langgraph"
set "MAMBA_ROOT_PREFIX=%ROOT%\micromamba"
set "KGENV=%ROOT%\envs\kali-hyp3"
set "TEMP=%ROOT%\opensource_work\tmp"
set "TMP=%ROOT%\opensource_work\tmp"
if not exist "%TEMP%" mkdir "%TEMP%"
if not exist "%MAMBA_ROOT_PREFIX%" mkdir "%MAMBA_ROOT_PREFIX%"

if exist "%MAMBA_ROOT_PREFIX%\micromamba.exe" goto :have_mamba
echo Downloading micromamba ...
curl.exe -L --fail -o "%MAMBA_ROOT_PREFIX%\micromamba.tar.bz2" https://micro.mamba.pm/api/micromamba/win-64/latest
if errorlevel 1 goto :fail
tar.exe -xjf "%MAMBA_ROOT_PREFIX%\micromamba.tar.bz2" -C "%MAMBA_ROOT_PREFIX%" Library/bin/micromamba.exe
if errorlevel 1 goto :fail
move /y "%MAMBA_ROOT_PREFIX%\Library\bin\micromamba.exe" "%MAMBA_ROOT_PREFIX%\micromamba.exe" >nul
rmdir /s /q "%MAMBA_ROOT_PREFIX%\Library"
del "%MAMBA_ROOT_PREFIX%\micromamba.tar.bz2"
:have_mamba

echo Creating the environment %KGENV% (about 500 MB to download) ...
"%MAMBA_ROOT_PREFIX%\micromamba.exe" create -y -p "%KGENV%" -f "%~dp0environment_windows.yml"
if errorlevel 1 goto :fail

call "%~dp0kg.cmd" check
echo.
echo Setup finished. Next:
echo   kg set-earthdata-token      (NASA Earthdata Login token, stored encrypted)
echo   kg set-cds-key              (only for the ERA5 variant)
echo   kg hyp3-search
exit /b 0

:fail
echo.
echo SETUP FAILED - see the messages above. Behind a university proxy see README_HYP3_WINDOWS.md, section "Network".
exit /b 1
