#!/bin/bash

function USAGE {
echo "Usage: "`basename $0`" <dir>";
echo ""
echo " Move IP4DI input and output files in <dir> for backup and post-"
echo " processing once the inversion has been run."
echo ""
}   # end function usage

#- if no input argument is given, display usage.
if [ $# -lt 1 ]; then
   USAGE
   exit 1
fi

#- create output directory
mkdir $1

#- mv output files
mv *.out $1
mv *.err $1
mv input*.mat $1
mv mesh*.mat $1
mv final*.mat $1
mv data*.dat $1
mv model*.dat $1
mv struct*.dat $1
mv ellipses.dat $1
mv info_vs_it*.txt $1

#- cp parameter files
cp -p forward_parameters.m $1
cp -p inversion_parameters.m $1
cp -p sensitivity_analysis.m $1

#- protect files from overwriting
chmod u-w $1/*
