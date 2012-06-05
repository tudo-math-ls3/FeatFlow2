#!/bin/bash

# Export output directory
export outputdir="fct_book/Shocktube_HartenEntropyFix"

IDISSIPATIONTYPE=( '1' '2' '3' )
for (( i=0;i<${#IDISSIPATIONTYPE[@]};i++)); do

# Export type of dissipation to be used
export idissipationtype=${IDISSIPATIONTYPE[${i}]}

NLEV=( '2' '3' '4' '5' '6' '7')
DSTEP=( '1.0e-3' '5.0e-4' '2.5e-4' '1.25e-4' '6.25e-5' '3.125e-5' )
for (( j=0;j<${#NLEV[@]};j++)); do

# Export level-dependent data
export nlev=${NLEV[${j}]}
export dstep=${DSTEP[${j}]}

# Low-order scheme
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_Loworder.dat) &

# Failsafe correction without main limiter
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_Failsafe-NoMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_Failsafe-WithMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_FailsafeAll-NoMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_FailsafeAll-WithMass.dat) &

# Wait until all background processes have finished
wait ${!}

LIMVAR1=( 'density' 'density,pressure' 'density,energy' )
LIMVAR=( 'r' 'rp' 'rE' )
for (( k=0;k<${#LIMVAR[@]};k++)); do

# Export settings for flux limiting
export limvar1=${LIMVAR1[${k}]}
export limvar=${LIMVAR[${k}]}

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT1-NoFailsafe-NoMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT1-WithFailsafe-NoMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT1-NoFailsafe-WithMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT1-WithFailsafe-WithMass.dat) &
done

# Wait until all background processes have finished
wait ${!}

LIMVAR1=( 'density' 'density' )
LIMVAR2=( 'pressure' 'energy' )
LIMVAR=( 'r,p' 'r,E' )
for (( k=0;k<${#LIMVAR[@]};k++)); do

# Export settings for flux limiting
export limvar1=${LIMVAR1[${k}]}
export limvar2=${LIMVAR2[${k}]}
export limvar=${LIMVAR[${k}]}

(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT2-NoFailsafe-NoMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT2-WithFailsafe-NoMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT2-NoFailsafe-WithMass.dat) &
(cd .. && ./flagship-pc64-nehalem-linux-intel-goto2 data/fct_book/1d/shocktube_LinFCT2-WithFailsafe-WithMass.dat) &
done

# Wait until all background processes have finished
wait ${!}

done

done
