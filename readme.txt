================================================================================
                                  FeatFlow 2.0
================================================================================

1. License And Copyright Information
2. Directory Structure
3. The Linux/Unix Build System
4. The Visual Studio + Intel Fortran Build System
   4.1 Third-Party Libraries
       4.1.1 Automatic Download And Unpack
       4.1.2 Manual Download And Unpack
       4.1.3 Compilation
   4.2 Compiling Applications
5. Final Remarks


================================================================================
1. License And Copyright Information
--------------------------------------------------------------------------------

Copyright (c) 1989-2014, Featflow Group
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
fully or in parts, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer. 
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


================================================================================
2. Directory Structure
--------------------------------------------------------------------------------

The FeatFlow 2 directory structure is organised as follows:

* applications
  > Contains the main applications like a poisson example, a heat equation
    example, cc2d for laminar Stokes/Navier-Stokes in 2D and others).

* bin
  > Basic shell scripts and makefile system

* docs
  > Contains scripts and input files for the documentation generator.

* kernel
  > The basic FeatFlow 2 kernel with all mathematical basics (finite element
    spaces, linear solvers, etc.)

* templates
  > Template scripts for different architectures for the creation of makefiles.

* tools
  > A set of auxiliary programs (e.g., creation of EPS files from meshes).

* tutorials
  > A set of tutorials demonstrating the usage of various components of the
    FeatFlow 2 kernel.


================================================================================
3. The Linux/Unix Build System
--------------------------------------------------------------------------------

To build a FeatFlow 2 application, you will need to perform the following steps:

1. Switch to the desired application directory (e.g., 'applications/poisson'):

      $ cd applications/poisson


2. Create a makefile script via a configure command. This will take a look at
   your system, detect your CPU and your compiler settings:

      $ ./configure

   Note: This is not the usual "configure" command which is provided by many
   Linux software packages. We wrote this configure on our own. It is a
   specialized shell/perl script and creates you a makefile "GNUmakefile".


3. Compile the application by entering

      $ make

   Note: During the first build, the build system will automatically download
   some basic libraries from the internet, untar/unzip the libraries and compile
   them (e.g., AMD, UMFPACK, LAPACK).


4. Run the desired application executable (e.g., 'poisson') by entering

      $ ./poisson


5. Optional: Start paraview and open one of the VTK files in the 'gmv'
   subdirectory. In the case of the 'poisson' applications, these files were
   represent the solution of a poisson equation on the unit square.


For further information on the build system, please refer to the 
'feat2buildsystem.pdf' inside the top-level directory.

================================================================================
4. The Visual Studio + Intel Fortran Build System
--------------------------------------------------------------------------------
To use FeatFlow 2 on Microsoft Windows, you will need:
 - Microsoft Visual Studio 2010 (or newer)
 - Intel Visual Fortran Compiler 12.0 (or newer)
 - Optional: Python 2.6 (or newer)


4.1 Third-Party Libraries
-------------------------
Before you can compile FeatFlow 2 or any of its applications, you will first
need to download and compile a set of third-party libraries, which are required
by the kernel modules.

There are two possibilies to download the required packages:
 - automatic download using a Python script
 - manual download and unpack


4.1.1 Automatic Download And Unpack
-----------------------------------
If you have an installation of the Python interpreter (version 2.6 or higher),
then you can simply execute the 'getlibs_win.py' script in the 'thirdparty'
directory - this will download and unpack all required third-party packages.
You need to be connected to the internet and allow internet access for the
Python interpreter in you firewall, of course.


4.1.2 Manual Download And Unpack
--------------------------------
If you do not have Python interpreter installed and do not plan to do so, you
can also download and unpack the packages manually. To do this, you need to
download all the download links which can be found in the 'links_win.txt' file
in the 'thirdparty' directory.

Moreover, you will need to unpack the archives manually. Unfortunately,
Microsoft Windows does not get shipped with a unpacker capable of unpacking
gzipped tarballs, so you also need to install an archive unpacking application.

If you do not already have an archive unpacker (e.g. WinZip, WinRAR or 7-Zip)
installed, we recommend the open-source unpacker 7-zip, which can be
downloaded at

                   http://www.7-zip.org/download

Once you have an unpacker ready to go, simply extract the contents of all
downloaded archives right into the 'thirdparty' directory.


4.1.3 Compilation
-----------------
To compile the third-party packages, you simply need to execute the
'buildlibs_win32.cmd' and/or 'buildlibs_win64.cmd', depending on whether you
plan to compile 32-bit and/or 64-bit versions of the applications, respectively.

Note:
If you are not using Visual Studio 2010 but a newer version, then these scripts
will fail. In this case, you need to open the 'thirdparty.sln' file from the
'thirdparty/visual_studio' directory in your desired Visual Studio version and
compile all required build configurations and platforms manually.

Once the compilation is finished, you are ready to compile FeatFlow 2 and its
applications.


4.2 Compiling Applications
--------------------------
Each of the applications in the 'applications' directory contains a 'solution'
file with the extension 'sln'. You can simply open this file in Visual Studio
and then compile the application binary for the desired build configuration and
platform.

There are 4 build configurations available:
- 'dbg'    : A standard (single-threaded) debug build
- 'dbg-omp': A debug build with OpenMP parallelisation
- 'opt'    : A (single-threaded) optimised build
- 'opt-omp': An optimised build with OpenMP parallelisation

Furthermore, two platforms are available:
- 'Win32': Builds a 32-bit Windows executable
- 'x64'  : Builds a 64-bit Windows executable

Note:
To compile builds for the 'x64' platform, you require the 64-bit versions of
both the Microsoft Visual C/C++ and the Intel Visual Fortran compiler.

Once the compilation is complete, you will find an executable in the application
directory, whose name is build up in the following way:

            <app-name>.if12-<platform>-<build-config>.exe

where:
 - <app-name> is the name of the application, e.g. 'poisson' or 'cc2d'
 - <platform> is either 'x86' (32-bit build) or 'x64' (64-bit build)
 - <build-config> is one of the four build configurations, e.g. 'dbg' or 'opt'


================================================================================
5. Final Remarks
--------------------------------------------------------------------------------

If you are using the software to compute results for a publication, we would
appreciate one of the following entries in the list of citations:

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


Good luck with your projects :)

  The Featflow-Team

Homepage:   http://www.featflow.de
EMail:      featflow@featflow.de
