#!/bin/bash

# This file is for searching for dependencies of USE-commands.
# Run it in with a list of Fortran .f90 files as argument.
# It will parse all .f90 files and search for USE commands.
# The results printed on screen can be copied to the makefile
# of the appropriate directory to tell MAKE which files depend
# on which files - so it's not anymore necessary to determine
# the order of the compilation of the .f90 files by hand.
# MAKE to recompile all dependend .f90 files!
#
# For this method to work properly, it's necessary that each
# module has the same name as the filename!

# creates dependency as: if file xxx.f90 is using module yyy then
# objdir/xxx.o moddir/xxx.mod: xxx.f90 moddir/yyy.mod 

for i in $@ 
  do 
  j=`basename $i .f90`
  iuse=`awk '/^[ ]*(U|u)(S|s)(E|e)/ {printf("$(MODDIR)/%s.mod\n",$2);}' $i | sort -u | tr -d "," | tr -s "\n" " " | tr -s "\'" " " `
  iinc=`awk '/^[ ]*(I|i)(N|n)(C|c)(L|l)(U|u)(D|d)(E|e)/ {print($2);}' $i | sort -u | tr -s "\n" " " | tr -s "\'" " " `
  printf "\$(OBJDIR)/${j}.o \$(MODDIR)/${j}.mod: ${i} ${iuse} ${iinc}\n"
done
