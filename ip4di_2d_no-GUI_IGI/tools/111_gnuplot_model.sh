#!/bin/bash
#
# Plot models stored in ASCII files in (x,z,value) format.
# NB: the 'reverse' option in yrange does not work anymore
#     with latest versions of Gnuplot (eg on Mac), use
#       plot * u 1:($z_max-\$2):3
#     instead (or execute the script from Linux).
#
# Francois Lavoue, Colorado School of Mines
# Last updated: March 2, 2016

file_in=$1.dat
cmplx_flag=$2
color_flag=$3

#title="$1"
#echo TITLE = $title

if [ $cmplx_flag -eq 3 ]; then
   suffix_out1=_amp

elif [ $cmplx_flag -eq 4 ]; then
   suffix_out1=_phase
fi

#write plot.gnu
echo "set term postscript eps enhanced color" > tmp_plot.gnu
echo "set output '$1$suffix_out1"_cs$color_flag".eps'" >> tmp_plot.gnu

#echo "set pm3d map" >> tmp_plot.gnu
#echo "set view map" >> tmp_plot.gnu
#echo "set pal gray" >> tmp_plot.gnu
#echo "set palette negative" >> tmp_plot.gnu   #reverse colorbar

# Matlab jet colormap
echo "set palette defined ( 0 '#000090',\\" >> tmp_plot.gnu
echo "                      1 '#000fff',\\" >> tmp_plot.gnu
echo "                      2 '#0090ff',\\" >> tmp_plot.gnu
echo "                      3 '#0fffee',\\" >> tmp_plot.gnu
echo "                      4 '#90ff70',\\" >> tmp_plot.gnu
echo "                      5 '#ffee00',\\" >> tmp_plot.gnu
echo "                      6 '#ff7000',\\" >> tmp_plot.gnu
echo "                      7 '#ee0000',\\" >> tmp_plot.gnu
echo "                      8 '#7f0000')  " >> tmp_plot.gnu

#echo "unset colorbox" >> tmp_plot.gnu

echo "set size ratio -1" >> tmp_plot.gnu

echo "set origin 0.0,0.0" >> tmp_plot.gnu
echo "set grid" >> tmp_plot.gnu
echo "show grid" >> tmp_plot.gnu

#echo "set title '"$title"'" >> tmp_plot.gnu

echo "set xlabel 'x (m)'" >> tmp_plot.gnu
echo "set ylabel 'z (m)'" >> tmp_plot.gnu

echo "set xrange [-0.025:1.225]" >> tmp_plot.gnu
echo "set yrange [-0.025:0.525] reverse" >> tmp_plot.gnu

##echo "set key bottom right" >> tmp_plot.gnu
##echo "set key width 3" >> tmp_plot.gnu     #noopaque ...

#echo "set parametric" >> tmp_plot.gnu
#echo "splot 'file_a.dat' using 1:2:3 with pm3d, 'curve.dat' u 1:2:(0.0) with lines"

if [ $cmplx_flag -eq 3 ]; then
# plot amplitude
   echo "set cblabel 'Resistivity, amplitude ({/Symbol W}.m)'" >> tmp_plot.gnu

   if [ $color_flag -eq 1 ]; then
   #same colormap as true model (1Hz)
      echo "set cbrange [275:375]" >> tmp_plot.gnu
      echo "set cbtics 275,25,375" >> tmp_plot.gnu
   fi

   echo "plot '"$PWD"/"$file_in"' u 1:2:(abs(\$3+{0.0,1.0}*\$4)) with image notitle" >> tmp_plot.gnu

elif [ $cmplx_flag -eq 4 ]; then
# plot phase
   echo "set cblabel 'Resistivity, phase (mrad)'" >> tmp_plot.gnu

   if [ $color_flag -eq 1 ]; then
   #same colormap as true model (1Hz)
      echo "set cbrange [0.4:1.8]" >> tmp_plot.gnu
      echo "set cbtics 0.5,0.25,1.75" >> tmp_plot.gnu
   fi

   echo "plot '"$PWD"/"$file_in"' u 1:2:(1000*atan(\$4/\$3)) with image notitle" >> tmp_plot.gnu

fi

# Gnuplot
echo "load '"$PWD"/tmp_plot.gnu'" | gnuplot -persist

# clean
rm tmp_plot.gnu
