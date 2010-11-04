# -*- mode: makefile -*-

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

# Intel-based Macs are quite tricky. Newer OSes claim to be 32bit but
# they are able to run 64bit code. Even worse, they generate 64bit
# code by default. However, if 32bit generation is forced then the
# linker must must be envoked with flag '-read_only_relocs suppress'.
ifeq ($(call match,$(ID),(pc|pc64)-.*-darwin-(g95|gcc)-.*),yes)
ifeq ($(strip $(BINARY)), 32)
LDFLAGS := -read_only_relocs suppress $(LDFLAGS)
endif
endif

# In this file, the third token has been set: operating system
# Set the flags accordingly.
TOKEN3 := 1

