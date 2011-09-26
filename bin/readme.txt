This directory contains auxiliary files for the compilation process.
The central file here is the shell-script GUESS_ID, which
creates an ID-string from the current computer architecture. This
ID-string is used in all build processes to uniquely name
subdirectories/executables for the current machine architecture,
so working on the same time on ALPHA as well as on INTEL as well
as on SUN workstations does not produce conflicts in the libraries/
executables - always the correct libraries/executables are used
on the current machine.
