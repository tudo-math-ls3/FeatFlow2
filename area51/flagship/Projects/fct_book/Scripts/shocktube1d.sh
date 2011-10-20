#!/bin/bash

#------------------------------------------------------------------------------
# Convergence study for Sod's shock tube problem in 1D
#
# dissipation: scalar dissipation proportional to the spectral radius of the Roe matrix
#              tensorial disspation (using Roe's linearization)
#              scalar dissipation of Rusanov type
#
# grid levels: 2-7
#
# time steps:  1e-3 - 3.125e-5
#
# discretizations: low-order, FCT-limiter with different limiting variables
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# Declaration of global parameters
#------------------------------------------------------------------------------

# Export output directory
export outputdir="fct_book/Shocktube1D_noEntropyfix"

# Export executable
export program="flagship-pc64-penryn-linux-intel-goto"

# Specify types of artificial dissipation
IDISSIPATIONTYPE=( '1' '2' '3' )

# Specify grid levels and time steps
NLEV=( '2' '3' '4' '5' '6' '7')
DSTEP=( '1.0e-3' '5.0e-4' '2.5e-4' '1.25e-4' '6.25e-5' '3.125e-5' )

#------------------------------------------------------------------------------
# Declaration of functions
#------------------------------------------------------------------------------

loworder()
{
    # Low-order scheme
    (eval $1/$program $1/data/fct_book/1d/shocktube_Loworder.dat) &
}

failsafe()
{
    # Failsafe correction for pressure and velocity without main limiter
    (eval $1/$program $1/data/fct_book/1d/shocktube_Failsafe-NoMass.dat) &
    (eval $1/$program $1/data/fct_book/1d/shocktube_Failsafe-WithMass.dat) &
}

failsafeall()
{
    # Failsafe correction for density, pressure and velocity without main limiter
    (eval $1/$program $1/data/fct_book/1d/shocktube_FailsafeAll-NoMass.dat) &
    (eval $1/$program $1/data/fct_book/1d/shocktube_FailsafeAll-WithMass.dat) &
}

linfct1_nofailsafe()
{
    LIMVAR1=( 'density' 'density,pressure' 'density,energy' )
    LIMVAR=( 'r' 'rp' 'rE' )
    for (( k=0;k<${#LIMVAR[@]};k++)); do
	
        # Export settings for flux limiting
	export limvar1=${LIMVAR1[${k}]}
	export limvar=${LIMVAR[${k}]}
	
	# Linearised FCT limiter with synchronisation and without failsafe feature
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT1-NoFailsafe-NoMass.dat) &
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT1-NoFailsafe-WithMass.dat) &
    done
}

linfct1_withfailsafe()
{
    LIMVAR1=( 'density' 'density,pressure' 'density,energy' )
    LIMVAR=( 'r' 'rp' 'rE' )
    for (( k=0;k<${#LIMVAR[@]};k++)); do
	
        # Export settings for flux limiting
	export limvar1=${LIMVAR1[${k}]}
	export limvar=${LIMVAR[${k}]}
	
	# Linearised FCT limiter with synchronisation and with
	# failsafe feature applied to pressure and velocity
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT1-WithFailsafe-NoMass.dat) &
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT1-WithFailsafe-WithMass.dat) &
    done
}

linfct1_withfailsafeall()
{
    LIMVAR1=( 'density' 'density,pressure' 'density,energy' )
    LIMVAR=( 'r' 'rp' 'rE' )
    for (( k=0;k<${#LIMVAR[@]};k++)); do
	
        # Export settings for flux limiting
	export limvar1=${LIMVAR1[${k}]}
	export limvar=${LIMVAR[${k}]}
	
	# Linearised FCT limiter with synchronisation and with
	# failsafe feature applied to density, pressure and velocity
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT1-WithFailsafeAll-NoMass.dat) &
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT1-WithFailsafeAll-WithMass.dat) &
    done
}

linfct2_nofailsafe()
{
    LIMVAR1=( 'density' 'density' )
    LIMVAR2=( 'pressure' 'energy' )
    LIMVAR=( 'r,p' 'r,E' )
    for (( k=0;k<${#LIMVAR[@]};k++)); do
	
        # Export settings for flux limiting
	export limvar1=${LIMVAR1[${k}]}
	export limvar2=${LIMVAR2[${k}]}
	export limvar=${LIMVAR[${k}]}
	
	# Linearised FCT limiter without synchronisation and without
	# failsafe feature
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT2-NoFailsafe-NoMass.dat) &
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT2-NoFailsafe-WithMass.dat) &
    done
}

linfct2_withfailsafe()
{
    LIMVAR1=( 'density' 'density' )
    LIMVAR2=( 'pressure' 'energy' )
    LIMVAR=( 'r,p' 'r,E' )
    for (( k=0;k<${#LIMVAR[@]};k++)); do
	
        # Export settings for flux limiting
	export limvar1=${LIMVAR1[${k}]}
	export limvar2=${LIMVAR2[${k}]}
	export limvar=${LIMVAR[${k}]}
	
        # Linearised FCT limiter without synchronisation and with
        # failsafe feature applied to pressure and velocity
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT2-WithFailsafe-NoMass.dat) &
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT2-WithFailsafe-WithMass.dat) &
    done
}

linfct2_withfailsafeall()
{
    LIMVAR1=( 'density' 'density' )
    LIMVAR2=( 'pressure' 'energy' )
    LIMVAR=( 'r,p' 'r,E' )
    for (( k=0;k<${#LIMVAR[@]};k++)); do
	
        # Export settings for flux limiting
	export limvar1=${LIMVAR1[${k}]}
	export limvar2=${LIMVAR2[${k}]}
	export limvar=${LIMVAR[${k}]}
	
        # Linearised FCT limiter without synchronisation and with
        # failsafe feature applied to density, pressure and velocity
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT2-WithFailsafeAll-NoMass.dat) &
	(eval $1/$program $1/data/fct_book/1d/shocktube_LinFCT2-WithFailsafeAll-WithMass.dat) &
    done
}

#------------------------------------------------------------------------------
# Main script
#------------------------------------------------------------------------------

# Loop over different dissipation types
for (( i=0;i<${#IDISSIPATIONTYPE[@]};i++)); do
    
    # Export type of dissipation to be used
    export idissipationtype=${IDISSIPATIONTYPE[${i}]}
    
    # Loop over different grid levels/time steps
    for (( j=0;j<${#NLEV[@]};j++)); do
	
        # Export level-dependent data
	export nlev=${NLEV[${j}]}
	export dstep=${DSTEP[${j}]}
	
        # Loop over all command line parameters
	for bench in "$@"; do
	    
	    case $bench in
		loworder)
		    echo "Computing low-order solution"
		    loworder "`pwd`"
		    ;;
		failsafe)
		    echo "Computing solution with failsafe as main limiter"
		    failsafe "`pwd`"
		    ;;
		failsafeall)
		    echo "Computing solution with failsafe for all variables as main limiter"
		    failsafeall "`pwd`"
		    ;;
		linfct1_nofailsafe)
		    echo "Computing solution with FCT1 without failsafe feature"
		    linfct1_nofailsafe "`pwd`"
		    ;;
		linfct1_withfailsafe)
		    echo "Computing solution with FCT1 with failsafe feature"
		    linfct1_withfailsafe "`pwd`"
		    ;;
		linfct1_withfailsafeall)
		    echo "Computing solution with FCT1 with failsafe feature for all variables"
		    linfct1_withfailsafeall "`pwd`"
		    ;;
		linfct2_nofailsafe)
		    echo "Computing solution with FCT2 without failsafe feature"
		    linfct2_nofailsafe "`pwd`"
		    ;;
		linfct2_withfailsafe)
		    echo "Computing solution with FCT2 with failsafe feature"
		    linfct2_withfailsafe "`pwd`"
		    ;;
		linfct2_withfailsafeall)
		    echo "Computing solution with FCT2 with failsafe feature for all variables"
		    linfct2_withfailsafeall "`pwd`"
		    ;;
		*)
		    echo "Invalid benchmark!"
		    ;;
	    esac
	done
        
        # Wait for unfinished jobs
        wait $!
    done
done
