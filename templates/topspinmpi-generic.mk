# -*- mode: makefile -*-

##############################################################################
# Environment flags for MPICH shipped with TOPSPIN MPI
# (see http://www.topspin.com/solutions/hpc.html)
##############################################################################

# If preprocessor switch -DENABLE_SERIAL_BUILD does not occur in compiler flags,
# a build for parallel execution is requested. MPI header and libraries possibly
# need to be added to INC, LIBS, LIBDIR.
ifeq (,$(findstring -DENABLE_SERIAL_BUILD ,$(APPONLYFLAGS) $(CFLAGSF90) ))

# In case MPI wrapper commands are used to compile code, no changes to INC, LIBS
# and LIBDIR are necessary. The MPI wrapper commands will take care of all the
# dirty stuff.
# If no MPI wrapper commands are to be used (sometimes they are not available at
# all or for different reasons one wishes to do so), add TopspinMPICH-related settings
# manually here.
ifeq ($(strip $(MPIWRAPPERS)), NO)

# Set MPICHHOME if unset
ifeq ($(strip $(MPICHHOME)),)
MPICHHOME = /usr/local/topspin/mpi/mpich
MESSAGE  := $(MESSAGE) \
	    echo '*** Warning: MPICHHOME unset. Has been set to $(MPICHHOME)'; \
	    echo '*** Warning: This has implications on include and library paths!';
endif

# Set include and library directory
MPIINC    = -I$(MPICHHOME)/include
MPILIBDIR = -L$(MPICHHOME)/lib64
MPILIBS   = -lmpich -lfmpich

endif  # MPIWRAPPERS=NO

endif  # MODE=PARALLEL aka -DENABLE_SERIAL_BUILD not set


# The settings needed to compile a FEAT2 application are "wildly" distributed
# over several files ((Makefile.inc and templates/*.mk) and if-branches 
# (in an attempt to reduce duplicate code and inconsistencies among all build 
# IDs that e.g. use the same MPI environment). Not having all settings at 
# *one* place entails the risk (especially in the event of setting up a new 
# build ID) that settings are incompletely defined. A simple typo in a matching 
# rule in Makefile.inc may prevent that the compiler and compiler command line
# flags are set. Compilation would fail with the most peculiar errors - if not
# the Makefile had been set up to catch such a case.
# Each build ID in FEAT2 has 6 tokens: architecture, cpu, operating system,
# compiler family, BLAS implementation, MPI environment. Whenever setting
# one of these, an according flag is set. They are named TOKEN1 up to TOKEN6.
# Before starting to actually compile a FEAT2 application, every Makefile
# generated by bin/configure checks whether all six tokens are set *for the 
# choosen build ID*. If not, the user gets an error message describing exactly 
# what information is missing, e.g. token 5 not set which means there is no
# information available which BLAS implementation to use and where to find the
# library.
#
# In this file, the sixth token has been set: MPI environment. 
# Set the flag accordingly.
TOKEN6 := 1
