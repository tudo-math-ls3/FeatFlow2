The following changes have been made in contrast to the origional
UMFPACK4 library:

-The "Source" folder has been renamed to "src".
-The "Include" folder has been renamed to "include".
-All unnecessary subfolders and the old makefile has been moved to
 the "orig" folder.
-In the umf_config.h file a new handling for BLAS routines was added:
 In case that the "UMF_WINDOWS" define exists, the BLAS routines are 
 expected to be in uppercase letters.
 (Internally this is checked the following way:
 Another include-file checks if the _WIN32 define of windows exists.
 If yes, it defines the "UMF_WINDOWS" define. If this exists,
 the BLAS function names are used in uppercase letters.)
-The "umf_multicompile_x.c" files have been added to support
 compiling without setting extra defines in the call of the compiler.
-The "umf4_f77wrapper.c" and "umf4_f77zrapper.c" have been moved
 to the "src" folder. All function names in this file have been
 converted to uppercase depending on the "UMF_WINDOWS" define,
 allowing the Windows compiler to link them to a Fortran program.
-The old makefile/GNUmakefile of the "Source" folder has been moved to
 "orig\Source".
-The makefile was rewritten in Feat-Style.

- make update will download the original UMFPACKv4.3.tar.gz from
  http://www.cise.ufl.edu/research/sparse/umfpack and perform all
  these changes. Our modified files are kept in the ./modsrc dir.