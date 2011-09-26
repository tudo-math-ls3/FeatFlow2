________________________________________________________________________
                    FEATFLOW 1.3 - libraries/
________________________________________________________________________

1. General structure
2. Libraries in the Microsoft Visual Studio IDE


                          1. General structure
------------------------------------------------------------------------
This directory contains the sourcecode for different libraries.
On one hand basic libraries can be found here like BLAS, LAPACK,...
On the other hand this directory also contains FEAT libraries
that are frequently used in all FEAT-applications - e.g. FEAT2D.

The source code of each library can be found in each of the
src/-subdirectories. The main directory contains a Makefile which
discribes how to compile the source files. For the FEAT-libraries
the sourcecode files are simply build with the compiler settings
in Globals.mk and then linked to a ".a"-library file.
This binary representations of the libraries are moved by the
makefiles into the directory "object/libraries" after compiling.

Remark: The "external libraries" in this source code package 
(LAPACK, BLAS, UMFPACK) are modified to be compiled in the same way.



          2. Libraries in the Microsoft Visual Studio IDE
------------------------------------------------------------------------
The library directories also contain projects for the 
Visual Studio 2003 IDE together with the Intel Fortran Compiler
for Windows. The Visual Studio project can be found in the
winxxxx/-subdirectory of those applications. The applications are
assembled from the source code in the following way:
 a) All source files are added to the project
 b) The project options are modified the following way:
      Check Array and String Bounds = No
      Calling Convention            = C, Reference
    The calling convention is changed in advance for future
    enhancements. Although this is not necessary for the current
    FeatFlow version, it's necessary if C-libraries (e.g. UMFPACK-4)
    is later added.
When compiling Windows libraries, the binary representation (.lib-
file) of the library is placed in the Debug/- or Release/ subdirectory
of the appropriate winxxxx/ folder of the library, as usual for
the Visual Studio IDE. The libraries are incorporated into the Windows
application using Project dependencies, so there's no need to collect
them in a separate directory like it's the case for Unix-like
environments.
