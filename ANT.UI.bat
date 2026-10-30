@echo off
rem Double-click to start ANT.UI. To get a desktop icon: right-click this file > Send to > Desktop (create shortcut).
cd /d "%~dp0"
where pythonw >nul 2>nul || (
    echo pythonw was not found on PATH. Install Python 3 from python.org and tick "Add python.exe to PATH".
    pause
    exit /b 1
)
start "" pythonw ANT.UI.py
