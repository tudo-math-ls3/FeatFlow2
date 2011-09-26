

           FEAT2 pre-compiled GotoBLAS libraries for Windows Readme
           --------------------------------------------------------

Index
-----
I.   Prologue
II.  Using the pre-compiled GotoBLAS libraries
III. Compiling your own GotoBLAS library
IV.  List of pre-compiled GotoBLAS


-------------------------------------------------------------------------------
I. Prologue
-------------------------------------------------------------------------------

If you want to compile FEAT2 under Windows using the Visual Studio project
files, then you may choose the GotoBLAS as combined BLAS/LAPACK library instead
of the reference netlib BLAS/LAPACK implementation.

Although the GotoBLAS implements only a small subset of the LAPACK routines,
neither the FEAT2 kernel nor its libraries (currently) use LAPACK routines
which are not part of the GotoBLAS implementation, so unless an application
needs additional LAPACK routines, it is possible to use the GotoBLAS under
Windows.

Unfortunately, the GotoBLAS library cannot be compiled under Visual Studio but
has to be compiled under cygwin, which can be quite time-consuming and even
impossible (for Win64) if you don't have the necessary compilers installed.

Therefore this directory contains a hand full of pre-compiled GotoBLAS
libraries for different architectures which can be used to compile FEAT2 using
the Visual Studio project files under Windows.

In this file you will find a list of all GotoBLAS libraries that have been
pre-compiled including some informations about the machine on which the library
has been compiled.


-------------------------------------------------------------------------------
II. Using the pre-compiled GotoBLAS libraries
-------------------------------------------------------------------------------
The FEAT2 Visual Studio project files are configured to automatically link
the 'winblaslapack.lib" (32 bit) or 'winblaslapack_x64.lib' (64 bit) in the
'libraries' directory.

If you want to use a specialised GotoBLAS library, then there are 2 steps you
have to make (assume that 'libgoto_XXXX.lib' is the GotoBLAS binary you want to
use):

1. Copy the 'libgoto_XXXX.lib' into the 'libraries' directory and rename it
   to 'winblaslapack.lib' (32 bit) or 'winblaslapack_x64.lib' (64 bit).

2. Copy the 'libgoto_XXXX.dll' into your Windows-directory (usually
   'C:\Windows').


-------------------------------------------------------------------------------
III. Compiling your own GotoBLAS library
-------------------------------------------------------------------------------
If you want to compile your own (32-bit) GotoBLAS library, you will need at 
least two things: a) cygwin and b) the GotoBLAS source code.

Cygwin is freeware and can be downloaded under

http://www.cygwin.com/


The GotoBLAS is open-source for academic use and can be downloaded under

http://www.tacc.utexas.edu/resources/software/#blas


The two basic steps that have to be done are:
1. Install cygwin including the 'devel' package
2. Unpack the GotoBLAS source and follow the instructions in the readme to
   compile it.


It becomes more interesting if you want to compile a 64-bit GotoBLAS library,
as the cygwin port of the GNU compilers does not offer x64 compilers :(
In this case (and this is currently the only possibility) you will also need
the PGI C-compiler (64-bit version, of course), which unfortunately is not
freeware. However, a free evaluation version of the necessary PGI C-compiler
can be downloaded under

http://www.pgroup.com/support/downloads.php

Once you have downloaded and installed the PGI compiler, you can select it in
the GotoBLAS Makefile and (try to) compile a 64-bit GotoBLAS library.


-------------------------------------------------------------------------------
IV. List of pre-compiled GotoBLAS
-------------------------------------------------------------------------------
This section contains a list of all pre-compiled GotoBLAS libraries which
can be found in this directory, including informations about the machine,
OS and CPU used to compile the library.

Remark:
To save some disk space, the libraries have been compressed into ZIP archives.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
= libgoto_barcelonap-r1.26
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GotoBLAS revision.........: 1.26
Instruction Set...........: x86
Compiled under............: Windows XP Professional SP3
Compiler used.............: PGI C-compiler v7.2-4 x86
Compiled on CPU...........: AMD Phenom X4 9650
SMP active................: Yes
Maximum SMP threads.......: 4

Additional Makefile.rule changes:
SMP = 1
MAX_THREADS = 4

Additional Notes:
This library should work on all AMD Phenom X4 processors under Win32.

-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
= libgoto_barcelonap_x64-r1.26
-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GotoBLAS revision.........: 1.26
Instruction Set...........: x64
Compiled under............: Windows Vista Business x64 SP1
Compiler used.............: PGI C-compiler v7.2-4 x64
Compiled on CPU...........: AMD Phenom X4 9650
SMP active................: Yes
Maximum SMP threads.......: 4

Additional Makefile.rule changes:
SMP = 1
MAX_THREADS = 4
BINARY64 = 1

Additional Notes:
This library should work on all AMD Phenom X4 processors under Win64.

-------------------------------------------------------------------------------
*End-Of-File*
