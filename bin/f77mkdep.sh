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
  iuse=`awk '/^[ ]*(U|u)(S|s)(E|e)/ {printf("$(MODDIR)/%s.mod\n",$2);}' $i | sort -u | tr -d "," | tr -s "\n" " " | tr -s "\'" " " `
  iinc=`awk '/^[ ]*(I|i)(N|n)(C|c)(L|l)(U|u)(D|d)(E|e)/ {print($2);}' $i | sort -u | tr -s "\n" " " | tr -s "\'" " " `
  printf "\$(OBJDIR)/${j}.o: ${i} ${iinc} ${iuse}\n"
done
