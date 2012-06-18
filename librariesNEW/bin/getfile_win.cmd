@echo off

rem Executed in the BIN/ subdirectory.
rem
rem The script downloads a .tar.gz / .tgz file and extracts it.

.\bin\wget %1 -O %~nx1
.\bin\gzip -f --stdout -d %~nx1 | .\bin\tar -xf -
