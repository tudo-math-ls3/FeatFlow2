#!/bin/bash

# Set path to GMV files
ucdpath=../flagship/out/zpinch_paper/

for ilim in Rusanov; do
    for itime in CN BE; do
	    
	echo "*********************************"
	echo "Discretization scheme: " $ilim
	echo "Time-stepping scheme:  " $itime
	echo "*********************************"

	for igmv in $(seq -f "%05g" 10 5 95); do

	    echo "*********************************"
            echo "Sequence number:       " $igmv
	    echo "*********************************"

	    for ilev in 1 2 3 4; do

		# Z-pinch in rz-coordinates
		zresCoarse=`echo 3*2^\($ilev-1\) | bc`
		rresCoarse=`echo 120*2^\($ilev-1\) | bc`
		
		zresFine=`echo 3*2^\($ilev\) | bc`
		rresFine=`echo 120*2^\($ilev\) | bc`
		
		ucdreffile=Rect\_$itime\_$ilim\_$zresCoarse"x"$rresCoarse/zpinch\_$zresCoarse"x"$rresCoarse\_rusanov.$igmv.gmv
		ucdfile=Rect\_$itime\_$ilim\_$zresFine"x"$rresFine/zpinch\_$zresFine"x"$rresFine\_rusanov.$igmv.gmv

		# Z-pinch in xy-coordinates
#		zresCoarse=`echo 60*2^\($ilev-1\) | bc`
#		rresCoarse=`echo 240*2^\($ilev-1\) | bc`
		
#		zresFine=`echo 60*2^\($ilev\) | bc`
#		rresFine=`echo 240*2^\($ilev\) | bc`
		
#		ucdreffile=FullCirc\_$itime\_$ilim\_$zresCoarse"x"$rresCoarse/zpinch\_cn\_$zresCoarse"x"$rresCoarse\_rusanov.$igmv.gmv
#		ucdfile=FullCirc\_$itime\_$ilim\_$zresFine"x"$rresFine/zpinch\_cn\_$zresFine"x"$rresFine\_rusanov.$igmv.gmv

#		echo $ucdpath$ucdreffile
#		echo $ucdpath$ucdfile

		./convstudy \
		    --ucdfile $ucdpath$ucdfile \
		    --ucdreffile $ucdpath$ucdreffile
		
		echo "---------------------------------"
	    done
	done
	
	echo "*********************************"
    done
done

for ilim in FCT; do
    for itime in CN BE; do
	    
	echo "*********************************"
	echo "Discretization scheme: " $ilim
	echo "Time-stepping scheme:  " $itime
	echo "*********************************"

	for igmv in $(seq -f "%05g" 10 5 95); do

	    echo "*********************************"
            echo "Sequence number:       " $igmv
	    echo "*********************************"

	    for ilev in 1 2 3 4; do

		# Z-pinch in rz-coordinates
		zresCoarse=`echo 3*2^\($ilev-1\) | bc`
		rresCoarse=`echo 120*2^\($ilev-1\) | bc`
		
		zresFine=`echo 3*2^\($ilev\) | bc`
		rresFine=`echo 120*2^\($ilev\) | bc`
		
		ucdreffile=Rect\_$itime\_$ilim\_$zresCoarse"x"$rresCoarse/zpinch\_$zresCoarse"x"$rresCoarse\_fct.$igmv.gmv
		ucdfile=Rect\_$itime\_$ilim\_$zresFine"x"$rresFine/zpinch\_$zresFine"x"$rresFine\_fct.$igmv.gmv

		# Z-pinch in xy-coordinates
#		zresCoarse=`echo 60*2^\($ilev-1\) | bc`
#		rresCoarse=`echo 240*2^\($ilev-1\) | bc`
		
#		zresFine=`echo 60*2^\($ilev\) | bc`
#		rresFine=`echo 240*2^\($ilev\) | bc`
		
#		ucdreffile=FullCirc\_$itime\_$ilim\_$zresCoarse"x"$rresCoarse/zpinch\_cn\_$zresCoarse"x"$rresCoarse\_rusanov.$igmv.gmv
#		ucdfile=FullCirc\_$itime\_$ilim\_$zresFine"x"$rresFine/zpinch\_cn\_$zresFine"x"$rresFine\_rusanov.$igmv.gmv

#		echo $ucdpath$ucdreffile
#		echo $ucdpath$ucdfile

		./convstudy \
		    --ucdfile $ucdpath$ucdfile \
		    --ucdreffile $ucdpath$ucdreffile
		
		echo "---------------------------------"
	    done
	done
	
	echo "*********************************"
    done
done
