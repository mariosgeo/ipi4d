#!/bin/bash

file_TI='demodata/sandbox_inversion_f1Hz/image_sandbox_251x601_h0.002_xzv.dat'
file_ellipses='ellipses.dat'

title="$1"
echo TITLE = $title

# choose to plot axes in m or wrt mesh indices
plot_indices=0

# plot topo?
plot_topo=0

# geometry
x0=-2.8240   # min(mesh.param_y)
y0=2.5       # min(mesh.param_x)

dx=5.6480
dy=5

# plot dimensions
xmin=-10.0
xmax=645.0
if [ $plot_topo == 1 ]; then
   ymin=-550  #25  #-450
   ymax=200

   dytic=50
   ytic_min=-550
   ytic_max=200
else
   ymin=-5.00
   ymax=72.50

   dytic=20
   ytic_min=0
   ytic_max=650
fi

# tics
dxtic=50
xtic_min=0
xtic_max=650

# color scale
cmin=-0.5
cmax=0.5

# aspect ratio
ratio=-1


if [ $plot_indices -eq 1 ]; then
# conversion factors from meters to mesh indices
# (beware of the extra parenthesis)
   dxtic=10
   dytic=5
   cx="-($x0))/$dx+1"
   cy="-($y0))/$dy+1"
   ratio=$(echo "($ymax-($ymin))/($xmax-($xmin))" | bc -l)
   xmin=$(echo "($xmin$cx" | bc -l)
   xmax=$(echo "($xmax$cx" | bc -l)
   ymin=$(echo "($ymin$cy" | bc -l)
   ymax=$(echo "($ymax$cy" | bc -l)
else   # plot axes in m
   cx=")"
   cy=")"
fi


# write plot.gnu
echo "set term postscript eps enhanced color" > tmp_plot.gnu
echo "set output 'ellipses.eps'" >> tmp_plot.gnu

echo "set palette gray negative" >> tmp_plot.gnu   #reverse colorbar

#echo "set style fill transparent solid 0.5 noborder" >> tmp_plot.gnu

echo "set size 1.4" >> tmp_plot.gnu
echo "set size ratio $ratio" >> tmp_plot.gnu

echo "set origin 0.0,0.0" >> tmp_plot.gnu
#echo "set grid" >> tmp_plot.gnu
#echo "show grid" >> tmp_plot.gnu

echo "set pointsize 0.10" >> tmp_plot.gnu

echo "set title '"$title"'" >> tmp_plot.gnu

echo "set x2label 'x (m)'" >> tmp_plot.gnu
echo "set ylabel 'z (m)'" >> tmp_plot.gnu

echo "set x2range [$xmin:$xmax]" >> tmp_plot.gnu
echo "set yrange  [$ymin:$ymax]" >> tmp_plot.gnu   #reverse
echo "set cbrange [$cmin:$cmax]" >> tmp_plot.gnu

echo "unset xtics" >> tmp_plot.gnu
echo "set x2tics $xtic_min,$dxtic,$xtic_max" >> tmp_plot.gnu
echo "set ytics  $ytic_min,$dytic,$ytic_max" >> tmp_plot.gnu

echo "unset colorbox" >> tmp_plot.gnu
##echo "set cbtics (1,2)" >> tmp_plot.gnu
##echo "set cblabel 'Sand type'" >> tmp_plot.gnu

# plot TI
if [ $plot_topo == 1 ]; then
   echo "plot '"$PWD"/"$file_TI"' u (\$1$cx:(-\$2$cy:3 w points pt 5 ps 0.5 lc palette notitle axes x2y1, \\" >> tmp_plot.gnu
else
   echo "plot '"$PWD"/"$file_TI"' u (\$1$cx:(\$2$cy:3 w image notitle axes x2y1, \\" >> tmp_plot.gnu
fi

# plot ellipses
if [ $plot_topo == 1 ]; then
   echo "      '"$PWD"/"$file_ellipses"' u (\$1$cx:(-\$2$cy w points pt 7 lw 0.10 lc 1 notitle axes x2y1" >> tmp_plot.gnu
else
   echo "      '"$PWD"/"$file_ellipses"' u (\$1$cx:(\$2$cy w points pt 7 lw 0.10 lc 1 notitle axes x2y1" >> tmp_plot.gnu
fi

# Gnuplot
echo "load '"$PWD"/tmp_plot.gnu'" | gnuplot -persist

# clean
rm tmp_plot.gnu
