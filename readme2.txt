________________________________________________________________________
                             FEATFLOW 2.0
________________________________________________________________________


1. Installation
2. Directory structure
3. General remarks
4. Makefile structure
5. Developing your own application
6. Final remarks


                          1. Installation
------------------------------------------------------------------------

For successful installation of Featflow in a Unix-like environment
you need the following software packages:

- F95-compatible compiler like "g95" or "gfortran"
- C-compatible compiler like "gcc"
- GNU-make ("make"/"gmake")
- a couple of standard Unix commands that belong to the
  coreutils package: env, cat, tee, sed, awk
- a sh-compatible shell interpreter like "bash"


To install Featflow, switch to the folder of the application you want
to compile, e.g.

 "cd application/poisson"
 
and start the configure script in this folder to create a makefile for
your system:

 "./configure"
 
This creates a file "GNUmakefile" that contains the parameters about
your compiler. In the next step, compile the program by typing

 "make"
  
If you encounter problems with the installation, you can check the
following:

* Please check that all files in the "./bin" subdirectory have the
  "executable" flag set. If necessary, perform a "chmod a+x ./bin/*".

* If configure fails to automatically retrieve a build id identifying
  your processor and compiler environment, you can try to fetch all
  possible build id's by typing
  
  "./configure --list-ids"
  
  The result is a list of possible build id's, e.g.
  
  "Valid build IDs (for pc64-opteron-linux hosts) match
  against the following regular expressions:
  * pc64-opteron-linux-.*-.*
  * pc64-opteron-linux-.*-blas.*
  * pc64-opteron-linux-.*-goto.*
  * pc64-opteron-linux-g95-.*
  * pc64-opteron-linux-gcc-.*
  * pc64-opteron-linux-intel-.*
  
  To choose a build id, just concatenate the strings and instruct
  configure to adapt to that build id. A possible configure with a 
  build id created from the above list is e.g.
  
  "./configure --id=pc64-opteron-linux-gcc-blas"
  
  which creates a makefile in this situation for a 64 bit opteron
  machine on a linux based system using the GCC compiler suite and a
  standard BLAS library.

* By running "make id" you should check if your machine is detected
  correctly and the compilers are found. If you have non-standard
  compiler (-settings) you can edit the Globals.mk file in the main
  Featflow directory to specify settings for youe system; see section
  "Makefile-structure" below for details. You can also use
  "make sysinfo" for a detailed information page about your
  computer system and the compiler settings that are chosen.
  
* If you interchanged files between Linux and Windows systems, check
  that the files (expecially the "feat.msg" file in each applications/-
  subdirectory and the .prm/.tri/.dat data files of each problem) are
  saved in Unix text format. Eventually you have to use "dos2unix" to
  convert these files.


For an overview about the other options valid for installation, type

  "./configure help"
  
for the configure script and

  "make help"

for the makefile of the application created by configure.


                      2. Directory structure
------------------------------------------------------------------------

The directory structure has changed in contrast to the original
manual due to the better handling of makefiles. The files and
applications are organized as follows:

./applications
->  Contains all Featflow-applications: cc2d, pp2d,...
    as well as the benchmark computation suite for the cylinder-flow
    benchmark (benchmark, bench_xxx)
   
./area51
->  Testing playground for applications that are still in development
    or for test versions.
   
./bin
->  Contains auxiliary shell scripts to support compiling on
    multiple computer architectures (Intel, Sun, Alpha,...)
    and the usage of include-files

./libraries
->  Contains the source code for all libraries - the general libraries
    (blas, lapack,...) as well as the FEAT libraries
    
./kernel
->  Contains the source of kernel modules that realise the core of the
    FEAT2 project.


                         3. General remarks
------------------------------------------------------------------------

a) It should be possible to compile this package on most computers with
   a Unix-like programming environment, e.g.
    - Sun Workstation
    - Alpha
    - Intel-compatible x86-computers with Linux
    - Windows-computers using the Cygwin Unix emulator package
      (http://www.cygwin.com/)
 
b) For the libraries as well as for the main applications, there are
   even project files provided for the Microsoft Visual Studio 2003 IDE 
   in combination with the Intel Fortran Compiler 8.0. We remark that 
   these files are experimental, and so we give no guarantee for them to
   work. The readme.txt-file in the applications/- and libraries/-
   directory describes a little bit more in detail these project files.












                        4. Makefile structure
------------------------------------------------------------------------

In this section we give a short overview about the structure of the
makefiles for that case, that compilation is not possible due to
problems with compilers/linkers and so on.

There exist 3 core-makefiles in the main directory of Featflow:
  Globals.mk    - compiler settings
  Rules_apps.mk - application-specific make-rules
  Rules_libs.mk - library-specific make-rules

These makefiles contain global preferences for all applications in the
Featflow package. The main entry point for the compilation of each
application is the file "Makefile" in each application/library
directory as well as in the Featflow installation directory itself.

All makefiles including the one in the top-level directory of Featflow
share the same structure:

* Define the relative or absolute path to the top-level directory of
  the Featflow installation with the $(FEATFLOW) environment variable.

* Include "Globals.mk" for compiler settings

* Define a list of source files of the application/library with the
  $(SRC) environment variable

* If the makefile belongs to an application, include "Rules_apps.mk".
  If the makefile belongs to an library, include "Rules_libs.mk".

Globals.mk : This file is the main entry point for the compiler
  settings. If anything goes wrong (e.g. Featflow cannot be compiled
  due to wrong compiler settings), check this file and include your
  compiler settings here.
  
Rules_apps.mk : This file contains the compiler rules for applications.
  Don't ever modify this file, since it's used by all makefiles
  belonging to applications in Featflow!
  
Rules_libs.mk : This file contains the compiler rules for libraries.
  Don't ever modify this file, since it's used by all makefiles
  belonging to libraries in Featflow!

The top-level "Makefile" in the Featflow installation directory controls
the installation itself and calls the Makefiles in the different
subdirectories depending on the type of installation. To compile an
individual application/library, it's also possible to switch into
the corresponding directory and type
  "make"
there. Make will automatically detect if anything is missing (e.g.
libraries) and compile it if necessary using the compiler options
defined in "Globals.mk".


                  5. Developing your own application
------------------------------------------------------------------------

To start working on your own application/problem,
* pick one of the provided applications (ccXd/ppXd/bouss) on which
  you want to base your developement and
* make a copy of it inside of the applications/-subdirectory.

If you place your copy outside of the applications/-subdirectory,
edit the Makefile of your application: Set the variable $(FEATFLOW) 
to point to the base directory of your FEATFLOW installation 
(absolute or relative). After this you can compile your application by 
the make command and use all the other make targets:

* make all     - compile the application, if needed library is missing
                 compile it too, it is the default target (same as
                 "make" without any target)
* make id      - to check the compiler settings
* make help    - short overwiev of these make targets
* make clean   - to remove the object files for current machine
* make purge   - to remove the object files and the executable for
                 current machine

The makefile will automatically use the needed libraries from the base
$(FEATFLOW) installation. It also allows you to keep
executables/object files for different machines simultaneously.

By editing the Makefile in your application you can specify the name
to be used for your executable by the variable $(EXEC).
The standart setting is
  EXEC=appname-$(ID)
which will create the executable with that name and $(ID) is replaced
by the machine identification string so for example the resulting
executable will be named appname-pc-opteron-linux if compiled on
machine wich is identified as pc-opteron-linux.


                           7. Final remarks
------------------------------------------------------------------------

Also check the other "readme.txt"-files in the application/libraries-
subdirectories. They explain a little bit more in detail about the
content, structure and organization of these directories.


Good luck with your projects :-)

  The Featflow-Team

Homepage:   http://www.featflow.de
EMail:      featflow@featflow.de


========================================================================

last but not least...

Copyright notice

Featflow Disclaimer
Copyright (c) 1989-2005, Featflow Group
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, fully or in parts, are permitted provided that the
following conditions are met:

* The Featflow distribution package contains source code of other
  software packages (LAPACK, BLAS, UMFPACK). These are distributed
  under their own licenses which have to be respected.

* Redistributions of source code must retain the above copyright
  notice, this list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.

If you are using the software to compute results for a publication, we
would appreciate one the following entry in the list of citations:

"FEATFLOW - Finite Element Software for the Incompressible
 Navier-Stokes Equations; http://www.featflow.de; featflow@featflow.de"

If you are using TeX, you can use the following BibTeX-entry:

@TechReport{ TurekBecker1999 ,
  author =     {Turek, S. and Becker, Ch.},
  title =      {{FEATFLOW} - {F}inite {E}lement {S}oftware for the
               {I}ncompressible {N}avier--{S}tokes {E}quations}},
 institution = {University of Dortmund},
 year =        {1999},
 type =        {User  Manual}, 
 address =     {\url{http://www.featflow.de/}, \url{featflow@featflow.de}},
}
