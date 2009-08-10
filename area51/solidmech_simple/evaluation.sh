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
