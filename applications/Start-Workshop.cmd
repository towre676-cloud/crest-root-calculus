@echo off
cd /d "%~dp0"
where py >nul 2>&1
if errorlevel 1 goto try_python
py -3 workshop.py
goto done
:try_python
where python >nul 2>&1
if errorlevel 1 goto missing
python workshop.py
goto done
:missing
echo Install Python 3.11 or newer from python.org, then open this file again.
:done
if errorlevel 1 pause
