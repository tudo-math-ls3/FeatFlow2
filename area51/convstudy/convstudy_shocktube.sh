#!/bin/bash

# Set path to GMV files
ucdpath=../flagship/out/shocktube_convstudy/quads/

for ifail in nofailsafe failsafe; do
    for ilim in Fct Rusanov; do
	for itime in BE CN; do
	    
	    echo "*********************************"
	    echo "Discretization scheme: " $ilim
	    echo "Time-stepping scheme:  " $itime
	    echo "Failsafe algorithm:    " $ifail
	    echo "*********************************"

	    for ilev in 1 2 3 4 5; do
		
		jlev=`expr $ilev + 1`
		
		#echo "Grid levels "$ilev" and "$jlev
		
		ucdreffile=shocktube\_$ifail\_$ilim\_$itime\_NLEV$ilev.gmv
		ucdfile=shocktube\_$ifail\_$ilim\_$itime\_NLEV$jlev.gmv
		
		./convstudy \
		    --ucdfile $ucdpath$ucdfile \
		    --ucdreffile $ucdpath$ucdreffile

		echo "---------------------------------"
	    done
	    
	    echo "*********************************"
	done
    done
done
