# -*- mode: makefile -*-

##############################################################################
# SUN Ultrasparc v8
#
##############################################################################

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-g95-blas.*),yes)
# invokes the GNU g95 compiler, local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/OS-sunos.mk
include $(TEMPLATESDIR)/COMP-g95.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
ifeq ($(call optimise), YES)
CFLAGSF77     := -mcpu=ultrasparc $(CFLAGSF77)
CFLAGSF90     := -mcpu=ultrasparc $(CFLAGSF90)
CFLAGSC       := -mcpu=ultrasparc $(CFLAGSC)
LDFLAGS       := -mcpu=ultrasparc $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-g95-perf.*),yes)
# invokes the GNU g95 compiler, uses LAM MPI, Perflib BLAS,
# local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/OS-sunos.mk
include $(TEMPLATESDIR)/COMP-g95.mk
include $(TEMPLATESDIR)/BLAS-sunperf.mk
CFLAGSC       := $(CFLAGSC) #-DNSUNPERF
ifeq ($(call optimise), YES)
CFLAGSF77     := -mcpu=ultrasparc $(CFLAGSF77)
CFLAGSF90     := -mcpu=ultrasparc $(CFLAGSF90)
CFLAGSC       := -mcpu=ultrasparc $(CFLAGSC)
LDFLAGS       := -mcpu=ultrasparc $(LDFLAGS)
endif
endif

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, uses LAM MPI, local BLAS,
# local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/OS-sunos.mk
include $(TEMPLATESDIR)/COMP-sunstudio.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
# Have 32bit binaries built on Sun
CFLAGSF77     := -xarch=native $(CFLAGSF77)
CFLAGSF90     := -xarch=native $(CFLAGSF90)
CFLAGSC       := -xarch=native $(CFLAGSC)
LDFLAGS       := -xarch=native $(LDFLAGS)
# -lmvec Needed? Performance benefit? I don't know (SB)
#
# If using -lmvec (some Vectorised Mathematical function) we also
# need socket and nsl libraries.
LIBS          := $(LIBS) -lmvec -lsocket -lnsl
# With static compiling, don't try to link against pthread library
# By the way, both dynamic and static compiling works! Just comment
# out the following three lines to get dynamic compiling again.
LDFLAGS       := $(LDFLAGS) -Bstatic
LIBS          := $(filter-out -lpthread, $(LIBS))
LIBS          := $(filter-out -ldl, $(LIBS))
MPILIBS       := $(filter-out -lpthread, $(MPILIBS))
MPILIBS       := $(filter-out -ldl, $(MPILIBS))
endif

ifeq ($(call match,$(ID),sun4u-sparcv8-sunos-sunstudio-perf.*),yes)
# Use Perflib BLAS,
include $(TEMPLATESDIR)/BLAS-sunperf.mk
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS)) -xlic_lib=sunperf
endif



##############################################################################
# SUN Ultrasparc v9
#
##############################################################################

ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/OS-sunos.mk
include $(TEMPLATESDIR)/COMP-sunstudio.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
# Have 64bit binaries built on Sun
CFLAGSF77     := $(CFLAGSF77)     -xarch=native64
CFLAGSF90     := $(CFLAGSF90)     -xarch=native64
CFLAGSC       := $(CFLAGSC)       -xarch=native64
LDFLAGS       := $(LDFLAGS)       -xarch=native64
# If using -lmvec (some Vectorised Mathematical function) we also
# need socket and nsl libraries.
LIBS          := $(LIBS) -lmvec -lsocket -lnsl
endif

ifeq ($(call match,$(ID),sun4u-sparcv9-sunos-sunstudio-perf.*),yes)
# Use Perflib BLAS,
include $(TEMPLATESDIR)/BLAS-sunperf.mk
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS)) -xlic_lib=sunperf
endif



##############################################################################
# SUN Whatever Solaris 10
#
##############################################################################

ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-sunstudio-.*),yes)
# invokes the Sun Studio Compiler suite, local umfpack
#
# not fully tested settings
MESSAGE := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: The settings for this build ID have not been thoroughly tested yet'; \
	    echo;
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/OS-sunos.mk
include $(TEMPLATESDIR)/COMP-sunstudio.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
# Have 64bit binaries built on Sun
CFLAGSF77     := $(CFLAGSF77)     -xarch=native64
CFLAGSF90     := $(CFLAGSF90)     -xarch=native64
CFLAGSC       := $(CFLAGSC)       -xarch=native64
LDFLAGS       := $(LDFLAGS)       -xarch=native64
# If using -lmvec (some Vectorised Mathematical function) we also
# need socket and nsl libraries.
LIBS          := $(LIBS) -lmvec -lsocket -lnsl
# Is it necessary to replace the way Sunperf Library is included?
# It works both ways! (SB)
LIBS          := $(filter-out -lsunperf, $(LIBS))
LIBS          := $(LIBS) -xlic_lib=sunperf
endif

ifeq ($(call match,$(ID),sun4v-sparcv9-sunos-gcc-.*),yes)
# invokes the GNU Compiler suite, local umfpack
#
# not fully tested settings
MESSAGE := $(MESSAGE) \
	    echo; \
	    echo '*** Warning: The settings for this build ID have not been thoroughly tested yet'; \
	    echo;
include $(TEMPLATESDIR)/sun4u-generic.mk
include $(TEMPLATESDIR)/OS-sunos.mk
include $(TEMPLATESDIR)/COMP-gcc.mk
# Umfpack install instructions advice to set this flag on Sun.
CFLAGSC       := -DNSUNPERF $(CFLAGSC)
ifeq ($(call optimise), YES)
#CFLAGSF77LIBS := -march=athlon64 $(CFLAGSF77LIBS)
#CFLAGSF77     := -march=athlon64 $(CFLAGSF77)
#CFLAGSF90     := -march=athlon64 $(CFLAGSF90)
#CFLAGSC       := -march=athlon64 $(CFLAGSC)
#CFLAGSCXX     := -march=athlon64 $(CFLAGSCXX)
#LDFLAGS       := -march=athlon64 $(LDFLAGS)
endif
endif


ifeq ($(call match,$(ID),pc64-opteronx2-linux-.*-goto2),yes)
# Use GotoBLAS 2
include $(TEMPLATESDIR)/BLAS-gotoblas2.mk
include $(TEMPLATESDIR)/OS-linux.mk
# Workaround for linker problem:
# ld: hidden symbol `__svml_cosf4' in libsvml.a(svml_stub_scos4.o) is referenced by DSO
ifeq ($(call match,$(ID),pc64-opteronx2-linux-intel-goto2.*),yes)
LIBS := $(LIBS) -lsvml
endif
endif

ifeq ($(call match,$(ID),pc64-opteron-linux-.*-.*-optmpich),yes)
# Use Allinea Opt MPI
include $(TEMPLATESDIR)/MPI-optmpich.mk
include $(TEMPLATESDIR)/OS-linux.mk
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
#
# Code compiles, but does not run with OPT 1.0rc4. Dynamic linker
# can not find libopt.so. (Not a FEAT2 problem, though.)
#
# Additionally, OPT 1.0rc4 needs additional arguments in MPI calls.
# (Not sure about OPT 1.0 and above!!) parallel.f90 and parallelsys.f90
# are not prepared for them any more. (SB, 2006-06-02)
#
# # # # # # # # # # #      WARNING !!!     # # # # # # # # # # #
CFLAGSF77     := $(CFLAGSF77) -assume 2underscore
CFLAGSF90     := $(filter-out -DMPICALLS_NO_UNDERSCORE, $(CFLAGSF90)) \
		-DMPICALLS_NO_UNDERSCORE -assume 2underscore
SRCEXTRA      := $(SRCEXTRA) mpir_iargc.f

# Previous settings
# CFLAGSF77 = -g -us -assume 2underscore
# CFLAGSF90 = $(CFLAGSF77) -module $(OBJDIR) -check bounds -traceback \
#             -DHAS_INTRINSIC_IARGC -DHAS_FLUSH -DMPICALLS_NO_UNDERSCORE
# CFLAGSC   = -g -fpstkchk
# LDFLAGS   =
# SRCEXTRA  := $(SRCEXTRA) optmpi_wrapper.c
# MPIINC    = -I/usr/local/opt/mpich/include
# MPILIBDIR = -L/usr/local/opt/mpich/lib
# MPILIBS   = -L/usr/local/opt/opt/lib -lopt -lfmpich -lmpich -lpthread -lrt
# #MPILIBS   = -lpmpich -lfmpich -lmpich -lpthread -lrt
endif


##############################################################################
# ARM Cortex-A9
#
##############################################################################

ifeq ($(call match,$(ID),armv7l-cortexa9-linux-gcc-.*),yes)
# invokes the GNU Compiler suite, local umfpack, local metis
include $(TEMPLATESDIR)/pc-generic.mk
include $(TEMPLATESDIR)/COMP-gcc.mk
ifeq ($(strip $(optimise)), YES)
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(CFLAGSF77LIBS)
CFLAGSF77     := -mcpu=cortex-a9 -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp \
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(CFLAGSF77)
CFLAGSF90     := -mcpu=cortex-a9 -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp \
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(CFLAGSF90)
CFLAGSC       := -mcpu=cortex-a9 -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp \
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(CFLAGSC)
CFLAGSCOPROC  := -mcpu=cortex-a9 -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp \
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(CFLAGSCOPROC)
CFLAGSCXX     := -mcpu=cortex-a9 -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp \
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(CFLAGSCXX)
LDFLAGS       := -mcpu=cortex-a9 -march=armv7-a -mtune=cortex-a9 -mfloat-abi=softfp \
                 -mfp=vfpv3-d16 -mfpu=vfp -fopenmp -funroll-loops \
                 -DENABLE_CORTEXA_OPTIMISATIONS $(LDFLAGS)
endif
endif








ifeq ($(call match,$(ID),ia64-itanium2-linux-intel-mkl-mpi),yes)
# invokes the intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/ia64-generic.mk
include $(TEMPLATESDIR)/OS-linux.mk
include $(TEMPLATESDIR)/COMP-intel64.mk
LIBS     := $(LIBS) -L/opt/MathKeisan/MKL/lib/64 -lguide -lmkl_lapack64 -lmkl
#LIBS     := $(LIBS) -L/nfs/home11/HLRS/xam/xamsbuij/intel/Compiler/11.0/074/mkl/lib/64 -lguide -lmkl_lapack -lmkl
SBB_LIBS := $(SBB_LIBS) -L/opt/MathKeisan/MKL/lib/64 -lguide -lmkl_lapack64 -lmkl
TOKEN5 := 1
TOKEN6 := 1
endif



##############################################################################
# Itanium 2 Dual Core
#
##############################################################################

ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-mkl.*),yes)
include $(TEMPLATESDIR)/OS-linux.mk
ifneq ($(strip $(INTEL_MKL_LIB)),)
# Split up string if multiple directories are given
# Note: Do not put whitespace between comma and the environment variable, because
#       if you do, a string like "-L /path/to/mkl" is the result and that string
#       won't make it into the command line.
LIBDIR   := $(LIBDIR) -L$(subst :, -L,$(INTEL_MKL_LIB))
endif
LIBS     := $(LIBS) -lmkl_lapack -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lguide -lpthread
SBB_LIBS := $(SBB_LIBS) -lmkl_lapack -lmkl_intel_thread -lmkl_core -lmkl_intel_lp64 -lguide -lpthread
TOKEN5 := 1
endif

ifeq ($(call match,$(ID),ia64-itanium2x2-linux-.*-.*-ompi),yes)
include $(TEMPLATESDIR)/OS-linux.mk
# invokes the Intel Fortran 77/90 compiler, local umfpack
include $(TEMPLATESDIR)/MPI-openmpi.mk
endif



##############################################################################
# NEC
#
##############################################################################

ifeq ($(call match,$(ID),sx6-none-superux-f90-keisan-mpi),yes)
TOKEN5 := 1
TOKEN6 := 1
include $(TEMPLATESDIR)/necsx-generic.mk
include $(TEMPLATESDIR)/COMP-necsx.mk
include $(TEMPLATESDIR)/OS-superux.mk
ifeq ($(call optimise), YES)
CFLAGSF77     := -sx6 $(CFLAGSF77)
CFLAGSF90     := -sx6 $(CFLAGSF90)
CFLAGSC       := -sx6 $(CFLAGSC)
LDFLAGS       := -sx6 $(LDFLAGS)
endif
INC       = 
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd
LIBS      = -lmetis -lumfpack -lamd -L/SX/opt/MathKeisan/lib -llapack -lblas
SBB_LIBS  = -L/SX/opt/MathKeisan/lib -llapack -lblas
MPILIBS   = -lmpi
endif

ifeq ($(call match,$(ID),sx8-none-superux-f90-.*-mpi),yes)
TOKEN6 := 1
include $(TEMPLATESDIR)/necsx-generic.mk
include $(TEMPLATESDIR)/COMP-necsx.mk
include $(TEMPLATESDIR)/OS-superux.mk
ifeq ($(call optimise), YES)
CFLAGSF77     := -sx8 $(CFLAGSF77)  #-ftrace
CFLAGSF90     := -sx8 $(CFLAGSF90) #-Wf~-L fmtlist map summary transform~ -ftrace
CFLAGSC       := -sx8 $(CFLAGSC) #-ftrace
LDFLAGS       := -sx8 $(LDFLAGS) #-ftrace
endif
INC       =
SBB_LIBS  = -llapack -lblas
endif


ifeq ($(call match,$(ID),sx8-none-superux-f90-keisan-mpi),yes)
TOKEN5 := 1
include $(TEMPLATESDIR)/OS-superux.mk
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd
#LIBDIR   := $(LIBDIR) -L/SX/opt/MathKeisan/lib
#LIBDIR   := $(LIBDIR) -L/SX/opt/mathkeisan/MK1_6/lib0/lib64
LIBS  := $(LIBS) -L/SX/opt/MathKeisan/lib -llapack -lblas
MPILIBS   = -lmpi
endif




ifeq ($(call match,$(ID),sx8-none-superux-f90-blas-mpi),yes)
TOKEN5 := 1
include $(TEMPLATESDIR)/OS-superux.mk
#BUILDLIB  = metis blas lapack umfpack amd
BUILDLIB  = metis umfpack amd blas lapack
#LIBDIR   := $(LIBDIR) -L/SX/usr/lib0/libp
MPILIBS   = -lmpi
MPILIBS   =
endif
