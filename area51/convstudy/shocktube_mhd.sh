#!/bin/bash

# Set path to GMV files
ucdpath=../flagship/out/shocktube_mhd/ConvergenceStudy/

for ilev in 1; do
		
    jlev=`expr $ilev + 1`
		
    ucdfile=shocktube\_mhd\_FCT\_CN\_NLEV$ilev.gmv
    ucdreffile=shocktube\_mhd\_FCT\_CN\_NLEV$jlev.gmv

    echo $ucdpath$ucdreffile
    echo $ucdpath$ucdfile

    ./convstudy \
	--ucdfile $ucdpath$ucdfile \
	--ucdreffile $ucdpath$ucdreffile
    
    echo "*********************************"
done
