@echo off

rem Executed in the BIN/ subdirectory.
rem
rem The script downloads a .tar.gz / .tgz file and extracts it.
rem
rem The script uses the Unix/Linux commands WGET, GZIP and TAR which are
rem the original binaries provided by the following Web-Page:
rem
rem     http://gnuwin32.sourceforge.net/


rem Check if the desired file already exists and, if so, skip the download

if not exist %~nx1 goto download

echo File '%~nx1' already exists - skipping download
goto unpack


:download

.\bin\wget %1 -O %~nx1


:unpack

rem Unpacking is currently disabled due to problems with file access rights
rem when using gzip to unpack the archives. Unpack the archives by hand.
rem echo Unpacking '%~nx1'...

rem .\bin\gzip -f --stdout -d %~nx1 | .\bin\tar -xf -
