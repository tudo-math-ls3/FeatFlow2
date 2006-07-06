#!/bin/bash

# This file is used to remove the dependencies created by
# f90mkdep.sh for those modules which correspond to libraries
# and, hence, are defined externally.

depfile=$1;

for i in $@
  do
  case ${i} in
  -I*)
    for file in `\ls ${i:2}/*.mod 2> /dev/null`;
      do
      modname=`basename ${file}`;
      sed -i -e "s^\$(MODDIR)/${modname}^${file}^g;" $depfile;
    done
    ;;
  *)
    ;;
  esac
done

