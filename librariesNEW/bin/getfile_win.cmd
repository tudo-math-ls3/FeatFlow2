@echo off

rem Executed in the BIN/ subdirectory.
rem
rem The script downloads a .tar.gz / .tgz file and extracts it.
rem
rem The script uses the Unix/Linux commands WGET, GZIP and TAR which are
rem the original binaries provided by the following Web-Page:
rem
rem     http://gnuwin32.sourceforge.net/

.\bin\wget %1 -O %~nx1
.\bin\gzip -f --stdout -d %~nx1 | .\bin\tar -xf -
