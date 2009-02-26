# -*- mode: sh -*-
#
############################################################################
#
# Purpose of this file:
#
# feat2cppswitches.h serves as the central file in FEAT 2.0 where all 
# preprocessor switches (or more properly named 'macro definitions'), 
# which are used throughout FEAT's kernel source files or libraries, 
# are extensively documented and, occassionally, activated by uncommenting
# them.
#
# By default all macro definitions should be inactive, i.e. the CVS version 
# of this file should not contain anything but shell-style comments and 
# empty lines.
# 
# The file will be read in by Makefile.inc and stripped for its comments.
# Any remaining non-empty strings will be turned into a C-style -D macro 
# definition and added to CFLAGSF90. Example: 
# A line containing the following string (without the leading hash character)
#    ENABLE_SERIAL_BUILD
# will be turned into 
#    CFLAGS := -DENABLE_SERIAL_BUILD $(CFLAGS)
# On compilation, any code enclosed in a condition like
#   #ifdef ENABLE_SERIAL_BUILD
#      print *, "hello world."
#   #endif
# will be activated. Conversely, any code in a negated condition like
#   #ifndef ENABLE_SERIAL_BUILD
#      print *, "goodbye world."
#   #endif
# will be stripped before passing the preprocessed source to the compiler.
#
# There is no need to add these C-style -D macro definitions to CFLAGSF77 
# as well as Fortran77 source files are not currently piped through the 
# preprocessor.
# It might be an option to add these macro definitions in future to 
# CFLAGSC and/or CFLAGSCXX, too. Additional C and C++ source files 
# as well libraries written in these languages would benefit from it.
#
# Changes to this file cause a recompilation of all Fortran 90 source 
# files (for this file being on the Makefile's prerequisite list of every
# object file). That might not be optimal as most macro definitions
# affect only few, if not merely a single source file. But it is the 
# cheapest way to implement it so far. The intention of this file being 
# mainly to document the various preprocessor switches that can be set in 
# FEAT, an occiasional complete recompilation of a FEAT application should
# be tolerable.
#
#
#
############################################################################
#
# Why does one want to use a preprocessor?
#
# The need for a preprocessor arises from the desire to have applications 
# run on multiple platforms. A preprocessor will ensure that platform 
# dependencies such as file I/O, interlanguage calling sequences, 
# conditional compilation and, last but not least, workarounds for 
# compiler bugs will be handled in a correct and consistent manner. By 
# providing easy to use directives, a preprocessor can minimise developer 
# burden.
#
# FEAT uses the C preprocessor and bypasses any possibly existing
# Fortran preprocessor. The reason is straight forward: Experience has
# taught us that availability of a Fortran preprocessor is by no means 
# guaranteed. In fact, the most advanced platforms often have barely a
# flawless compiler full supporting the Fortran90/95 standard. A C 
# preprocessor, be it the GNU C preprocessor, is usually readily 
# available, though.
#



############################################################################
# General settings
