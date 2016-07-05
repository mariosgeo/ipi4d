#!/bin/bash

il_Cinit=6   #6 if IGI, 4 otherwise

LAGRN_LIST="1e-05 0.0001 0.001 0.01 0.1 1 10"    #1e-05

xoff1=0.01
yoff1=0.00002

xoff2=-0.01
yoff2=-0.0001

base_file_in=info_vs_it_lagrn

#write tmp_plot.gnu
echo "set term postscript eps enhanced color" > tmp_plot.gnu
echo "set output 'Lcurve.eps'" >> tmp_plot.gnu

echo "set xzeroaxis linetype -1 linewidth 0.5 " >> tmp_plot.gnu
echo "set yzeroaxis linetype -1 linewidth 0.5 " >> tmp_plot.gnu

#echo "set size ratio -1" >> tmp_plot.gnu
echo "set size 0.9,0.7" >> tmp_plot.gnu

echo "set origin 0.0,0.0" >> tmp_plot.gnu
echo "set grid" >> tmp_plot.gnu
echo "show grid" >> tmp_plot.gnu

echo "set pointsize 1.5" >> tmp_plot.gnu

echo "set title '"$title"'" >> tmp_plot.gnu

echo "set xlabel 'C_M/C_{init}'" >> tmp_plot.gnu
echo "set ylabel 'C_D/C_{init}'" >> tmp_plot.gnu

#echo "set log x" >> tmp_plot.gnu
#echo "set log y" >> tmp_plot.gnu

echo "set xrange [ 0:0.5]" >> tmp_plot.gnu
echo "set yrange [-0.0005:0.003]" >> tmp_plot.gnu

##echo "set xtics 1,1,$nsteps" >> tmp_plot.gnu
##echo "set ytics (0.01,0.025,0.05,0.1,0.25,0.5,1.0)" >> tmp_plot.gnu

##echo "set key bottom right" >> tmp_plot.gnu
##echo "set key width 3" >> tmp_plot.gnu     #noopaque ...


# init. data file
rm tmp_in.txt
il=1

# loop over Lagrangian values
for lagrn in $LAGRN_LIST
do

 file_in="$base_file_in"$lagrn.txt

 # read file_in
 nl=$(wc <$file_in | awk '{print $1}')
 lagrn2=$(awk <$file_in "NR==$nl {print \$7}")
 nit=$(awk <$file_in "NR==$nl {print \$1}")
 C_D=$(awk <$file_in "NR==$nl {print \$3}")
 C_M=$(awk <$file_in "NR==$nl {print \$4}")
 C_init=$(awk <$file_in "NR==$il_Cinit {print \$2}")

 # reminder: cols of tmp_in
 # 1: iteration nb
 # 2: C_tot(m_k)
 # 3: C_D(m_k)
 # 4: C_M(m_k)
 # 5: lagrn*C_M(m_k)
 # 6: lagrn*C_M(m_k)/C_D(m_k)
 # 7: lagrn
 # 8: normalized RMS

 # checks
 echo " "
 echo "========================"
 echo "lagrn = $lagrn, $lagrn2"
 echo "N_it = $nit"
 echo "C_D = $C_D"
 echo "C_M = $C_M"

 # make data file
 echo "$C_M $C_D" >> tmp_in.txt

 # set labels (lagrn values and nit)
 #echo "set label $il '{/Symbol l}=$lagrn, n_{it}=$nit' at first 1.1*$C_M,1.1*$C_D/$C_init font 'Roman,10'" >> tmp_plot.gnu   # LOG

 if [ $il -eq 2 ] || [ $il -eq 12 ] ; then
    xl=$(echo "1.0*($C_M/$C_init+$xoff2)" | bc -l)
    yl=$(echo "1.0*($C_D/$C_init+$yoff2)" | bc -l)
    echo "set label $il '{/Symbol l}=$lagrn, n_{it}=$nit' at first $xl,$yl font 'Roman,10'" >> tmp_plot.gnu
 
 else
    xl=$(echo "1.0*($C_M/$C_init+$xoff1)" | bc -l)
    yl=$(echo "1.0*($C_D/$C_init+$yoff1)" | bc -l)
    echo "set label $il '{/Symbol l}=$lagrn, n_{it}=$nit' at first $xl,$yl font 'Roman,10'" >> tmp_plot.gnu
 fi
 il=$(($il+1))

done   #end loop over lagrn

# plot
echo "plot '"$PWD"/tmp_in.txt' u (\$1/$C_init):(\$2/$C_init) with linespoints lw 2 lc 8 notitle" >> tmp_plot.gnu

# Gnuplot
echo "load '"$PWD"/tmp_plot.gnu'" | gnuplot -persist

# clean
#rm tmp*
