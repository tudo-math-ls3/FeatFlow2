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


function calc_err_red () {
  local pattern;
  local file;
  local old;
  local i;
  local red;
  pattern=$1;
  file=$2;

  if test -z "$pattern" -o "$pattern" = "-h" -o "$pattern" = "--help"; then
    echo $progname": Script which calculates reducion rate.";
    echo;
    echo "Usage:";
    echo "   $progname <pattern to match> <infile>";
    echo;
    echo "e.g.: $progname \"^abs. L2 error:\" results001.log";
    exit;
  fi

  old="";
#  echo "# * matches for <$pattern>:"
  for i in $(grep "${pattern}" "${file}" | awk '{print $(NF)}');
  do
    red="-";
    if test "$old" != ""; then
      red=$(printf "%6.3f" $(math ${old} / ${i}));
    fi
    printf "%8.4e   %5s\n" ${i} ${red}
    old=${i};
  done
}

if test $# -ne 1 -o "$1" = "-h" -o "$1" = "--help"; then
    usage;
    exit
fi

echo "--------------------------------"
echo "L2 error:"
grep "L2 error for  u:" $* | awk '{print $(NF)}'
echo "L2 error reduction:"
calc_err_red "L2 error for  u:" $*
echo "--------------------------------"
echo "H1 error:"
grep "H1 error for  u:" $* | awk '{print $(NF)}'
echo "H1 error reduction:"
calc_err_red "H1 error for  u:" $*
echo "--------------------------------"

echo "#DOF:"
grep "Number of DOF:" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Iterations:"
grep "Number of iterations" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Convergence rate:"
grep "Convergence rate" $* | awk '{print $(NF)}'
echo "--------------------------------"

echo "Convergence rate (asymptotic):"
grep "Asymptotic convergence rate" $* | awk '{print $(NF)}'
echo "--------------------------------"

diverged=$(grep "Solution diverging!" $*)
if test "$diverged" != ""; then
  echo $diverged
  echo "--------------------------------"
fi

#echo "U1 Displacement"
#grep "U1(X,Y):" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "U2 Displacement"
#grep "U2(X,Y):" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Abs error"
#grep "Abs-error for U:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Strain eps11 "
#grep "eps11:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Strain eps22 "
#grep "eps22:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Strain eps12 "
#grep "eps12:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Stress sigma11 "
#grep "sigma11:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Stress sigma22 "
#grep "sigma22:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Stress sigma33 "
#grep "sigma33:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Stress sigma12 "
#grep "sigma12:" $* | awk '{print $(NF)}'
#echo "--------------------------------"
#
#echo "Norm of Devsigma "
#grep "|devsigma|:" $* | awk '{print $(NF)}'
#echo "--------------------------------"