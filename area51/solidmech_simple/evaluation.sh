#!/bin/bash
#
# Usage:
#   evaluation <infile>
#
# Author: sven.buijssen@math.uni-dortmund.de,hilmar.wobker@math.uni-dortmund.de  2007/02/09
#

progname=$0

function usage ()
{
    echo $progname": Script which scans a special kind of log file"
    i=1; while test $i -le ${#progname}; do echo -n " "; i=$(expr $i + 1); done
    echo "  and outputs the crucial information."
    echo
    echo "Usage:"
    echo "   $progname <infile>"
    echo
    echo "e.g.: $progname results001.log"
    echo
}


function math ()
{
    export awktmp="$*";
    echo "" | awk "{ print $awktmp }";
    unset awktmp
}


if test $# -ne 1 -o "$1" = "-h" -o "$1" = "--help"; then
    usage;
    exit
fi

echo "--------------------------------"
echo "L2 error"
grep "L2-error for U:" $* | awk '{print $(NF)}'

echo "--------------------------------"
echo "H1 error"
grep "H1-error for U:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "#DOF"
grep "Number of DOF:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Iterations"
grep "Iterations " $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Rate of Convergence"
grep "Rate of Convergence" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Rate of Asymptotic Convergence"
grep "Rate of Asymptotic Convergence:" $* | awk '{print $(NF)}'
echo "--------------------------------"

grep "Solution diverging!" $*
echo "--------------------------------"

echo "U1 Displacement"
grep "U1(X,Y):" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "U2 Displacement"
grep "U2(X,Y):" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Abs error"
grep "Abs-error for U:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Strain eps11 "
grep "eps11:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Strain eps22 "
grep "eps22:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Strain eps12 "
grep "eps12:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Stress sigma11 "
grep "sigma11:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Stress sigma22 "
grep "sigma22:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Stress sigma33 "
grep "sigma33:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Stress sigma12 "
grep "sigma12:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Norm of Devsigma "
grep "|devsigma|:" $* | awk '{print $(NF)}'
echo "--------------------------------"