#!/bin/bash

#file_TI=sandbox_structure_v2_251x601_h0.002.bin
file_TI=mgradient_init_sandbox_v2_251x601_h0.002_xzv.dat
file_ellipses=ellipses.dat

title="$1"
echo TITLE = $title

#write plot.gnu
#echo "set term postscript eps enhanced color" > tmp_plot.gnu
#echo "set output 'ellipses $title.eps'" >> tmp_plot.gnu

#echo "set xzeroaxis linetype -1 linewidth 0.5 " >> tmp_plot.gnu
#echo "set yzeroaxis linetype -1 linewidth 0.5 " >> tmp_plot.gnu

#echo "set pm3d map" >> tmp_plot.gnu
#echo "set view map" >> tmp_plot.gnu
echo "set pal gray" >> tmp_plot.gnu
echo "set palette negative" >> tmp_plot.gnu   #reverse colorbar

echo "set size ratio -1" >> tmp_plot.gnu

echo "set origin 0.0,0.0" >> tmp_plot.gnu
echo "set grid" >> tmp_plot.gnu
echo "show grid" >> tmp_plot.gnu

echo "set pointsize 0.25" >> tmp_plot.gnu

echo "set title '"$title"'" >> tmp_plot.gnu

echo "set xlabel 'x (m)'" >> tmp_plot.gnu
echo "set ylabel 'z (m)'" >> tmp_plot.gnu

echo "set xrange [-0.025:1.225]" >> tmp_plot.gnu
echo "set yrange [-0.025:0.525] reverse" >> tmp_plot.gnu

echo "set cbrange [-1.5e7:1.5e7]" >> tmp_plot.gnu
#echo "set cbtics (1,2)" >> tmp_plot.gnu
#echo "set cblabel 'Sand type'" >> tmp_plot.gnu

echo "unset colorbox" >> tmp_plot.gnu


##echo "set xtics 1,1,$nsteps" >> tmp_plot.gnu
##echo "set ytics (0.01,0.025,0.05,0.1,0.25,0.5,1.0)" >> tmp_plot.gnu

##echo "set key bottom right" >> tmp_plot.gnu
##echo "set key width 3" >> tmp_plot.gnu     #noopaque ...

#echo "set parametric" >> tmp_plot.gnu
#echo "splot 'file_a.dat' using 1:2:3 with pm3d, 'curve.dat' u 1:2:(0.0) with lines"

# plot TI
echo "plot '"$PWD"/"$file_TI"' u 1:2:(1+\$3) with image palette notitle, \\" >> tmp_plot.gnu

# plot ellipses
echo "      '"$PWD"/"$file_ellipses"' u 1:2 with points lw 0.5 pt 7 lc 1 notitle" >> tmp_plot.gnu


echo "load '"$PWD"/tmp_plot.gnu'" | gnuplot -persist

# clean
rm tmp_plot.gnu
