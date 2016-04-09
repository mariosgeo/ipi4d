#!/bin/bash

base_file_in=model_final_lagrn

# loop over Lagrangian values
for lagrn in 0.0001 0.001 0.01 0.1 1 10    #1e-05
do

  file_in="$base_file_in"$lagrn

  for color_flag in 1 2
  do
  # 1 -> plot with same color scale as true model
  # 2 -> automatic color scale

    for cmplx_flag in 3 4
    do
    # 3->amplitude
    # 4->phase
      ../111_gnuplot_model.sh $file_in $cmplx_flag $color_flag
    done

  done   #end color_flag

done   #end loop over lagrn

# clean
rm tmp*
