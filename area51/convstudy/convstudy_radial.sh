#!/bin/bash

# Set path to GMV files
ucdpath=../flagship/out/radial_convstudy/quads/

for ifail in nofailsafe failsafe; do
    for ilim in Fct Rusanov; do
	for itime in BE CN; do
	    
	    echo "*********************************"
	    echo "Discretization scheme: " $ilim
	    echo "Time-stepping scheme:  " $itime
	    echo "Failsafe algorithm:    " $ifail
	    echo "*********************************"

	    for ilev in 4 5 6 7 8; do
		
		jlev=`expr $ilev + 1`
		
		#echo "Grid levels "$ilev" and "$jlev
		
		ucdreffile=radial\_$ifail\_$ilim\_$itime\_NLEV$ilev.gmv
		ucdfile=radial\_$ifail\_$ilim\_$itime\_NLEV$jlev.gmv
		
		./convstudy \
		    --ucdfile $ucdpath$ucdfile \
		    --ucdreffile $ucdpath$ucdreffile

		echo "---------------------------------"
	    done
	    
	    echo "*********************************"
	done
    done
done
