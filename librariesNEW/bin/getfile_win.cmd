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

echo Unpacking '%~nx1'...
.\bin\7za x -so %~nx1 2>&3 | .\bin\7za x -si -ttar -aos > nul