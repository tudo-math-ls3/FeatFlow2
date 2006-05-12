#!/bin/bash

# This file is for searching for dependencies of INCLUDE-files.
# Run it in with a list of Fortran .f files as argument.
# It will parse all .f files and search for INCLUDE commands.
# The results printed on screen can be copied to the makefile
# of the appropriate directory to tell MAKE which files depend
# on which include files - so editing an include file will force
# MAKE to recompile all dependend .f files!

for i in $@ 
  do 
  j=`basename $i .f`
  printf "\$(OBJDIR)/%s.o: %s.f" $j $j
  awk '/^[ ]*(I|i)(N|n)(C|c)(L|l)(U|u)(D|d)(E|e)/ {printf(" %s",$2);}' $i | tr -d \' 
  printf "\n" 
done
