#!/bin/bash
#
# Program that moves IP4DI output files in a directory.
# To use it, type in your Terminal
# ./pack_output_IP4DI.sh directory_name

# create output directory
mkdir $1

# mv output files
mv log.out $1
mv debug.out $1
mv mesh*.mat $1
mv data*.dat $1
mv model*.dat $1
mv struct*.dat $1
mv ellipses.dat $1
mv info_vs_it*.txt $1

# cp parameter files
cp -p forward_parameters.m $1
cp -p inversion_parameters.m $1
cp -p sensitivity_analysis.m $1
