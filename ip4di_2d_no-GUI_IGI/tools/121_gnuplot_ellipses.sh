#!/bin/bash

file_TI='demodata/sandbox_inversion_f1Hz/image_sandbox_251x601_h0.002_xzv.dat'
file_ellipses=ellipses.dat

title="$1"
echo TITLE = $title

#write plot.gnu
echo "set term postscript eps enhanced color" > tmp_plot.gnu
echo "set output 'ellipses.eps'" >> tmp_plot.gnu

echo "set palette gray negative" >> tmp_plot.gnu   #reverse colorbar

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

echo "unset colorbox" >> tmp_plot.gnu
##echo "set cbtics (1,2)" >> tmp_plot.gnu
##echo "set cblabel 'Sand type'" >> tmp_plot.gnu

# plot TI
echo "plot '"$PWD"/"$file_TI"' u 1:2:3 w image notitle, \\" >> tmp_plot.gnu

# plot ellipses
echo "      '"$PWD"/"$file_ellipses"' u 1:2 with points lw 0.5 pt 7 lc 1 notitle" >> tmp_plot.gnu

# Gnuplot
echo "load '"$PWD"/tmp_plot.gnu'" | gnuplot -persist

# clean
rm tmp_plot.gnu
