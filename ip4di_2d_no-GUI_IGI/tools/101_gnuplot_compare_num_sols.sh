#!/bin/bash

freq=1
input=$1

dir1=data-init_sandbox_mesh-inv
dir2=data-init_sandbox_mesh-inv_f1Hz

file1=forward_sandbox_mesh-inv_data-init_f"$freq"Hz.dat
file2=forward_sandbox_mesh-inv_data-init_f"$freq"Hz.dat

title1='Analytical init. data'
title2='Numerical init. data'
title3='analytical vs numerical init. data'

file_out=forward_sandbox_data-init_ana-vs-num_f"$freq"Hz

# recast both datasets in same file for comparison
awk <$dir1/$file1 '{print $9,$10}' >tmp1
awk <$dir2/$file2 '{print $9,$10}' >tmp2

paste tmp1 tmp2 >tmp.dat

#write plot.gnu
echo "set term postscript eps enhanced color" > tmp_plot.gnu

echo "set xzeroaxis linetype -1 linewidth 0.5 " >> tmp_plot.gnu
echo "set yzeroaxis linetype -1 linewidth 0.5 " >> tmp_plot.gnu

echo "set size 1,0.8" >> tmp_plot.gnu

echo "set origin 0.0,0.0" >> tmp_plot.gnu
echo "set grid" >> tmp_plot.gnu
echo "show grid" >> tmp_plot.gnu

echo "set pointsize 1" >> tmp_plot.gnu

echo "set ylabel 'Apparent resistivity (Ohm.m)'" >> tmp_plot.gnu
#echo "set yrange [460:500]" >> tmp_plot.gnu

## plot wrt data nb
echo "set xlabel 'Data nb'" >> tmp_plot.gnu

#if [ $input==1 ]; then
  # plot real parts

#  echo "set term x11 0" >> tmp_plot.gnu

  echo "set output '"$file_out"_real.eps'" >> tmp_plot.gnu

  echo "plot '"$PWD"/"tmp.dat"' u 1 w lines lw 2 lt 1 lc 1 title '"$title1" (real part)', \\" >> tmp_plot.gnu
  echo "     '"$PWD"/"tmp.dat"' u 3 w lines lw 2 lt 2 lc 2 title '"$title2" (real part)'    " >> tmp_plot.gnu

#elif [ $input==2 ]; then
  # plot imag. parts

#  echo "set term x11 1" >> tmp_plot.gnu

  echo "set output '"$file_out"_imag.eps'" >> tmp_plot.gnu

  echo "plot '"$PWD"/"tmp.dat"' u 2 w lines lw 2 lt 1 lc 1 title '"$title1" (imag. part)', \\" >> tmp_plot.gnu
  echo "     '"$PWD"/"tmp.dat"' u 4 w lines lw 2 lt 2 lc 2 title '"$title2" (imag. part)'    " >> tmp_plot.gnu

#elif [ $input==3 ]; then
  # plot relative difference

#  echo "set term x11 2" >> tmp_plot.gnu

  echo "set output '"$file_out"_diff.eps'" >> tmp_plot.gnu

  echo "i = {0.0,1.0}" >> tmp_plot.gnu
  echo "set log y" >> tmp_plot.gnu
#  echo "set yrange [-1e-8:1e-8]" >> tmp_plot.gnu
  echo "set ylabel 'Relative difference |1-2|/|1|'" >> tmp_plot.gnu
  echo "plot '"$PWD"/"tmp.dat"' u (abs((\$1+i*\$2)-(\$3+i*\$4))/abs(\$1+i*\$2)) w lines lw 2 lt 1 lc 1 title 'Relative difference "$title3"'" >> tmp_plot.gnu

#fi

echo "load '"$PWD"/"tmp_plot.gnu"'" | gnuplot -persist

# clean
rm tmp*
