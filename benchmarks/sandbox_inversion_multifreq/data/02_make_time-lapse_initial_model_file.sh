#!/bin/bash

ls_file_in=ls_model_files.txt
file_out=model_sandbox_init_multifreq.mat

nfiles=$(wc <$ls_file_in | awk '{print $1}')
echo $nfiles > tmp_nfiles.txt

for ((ifile=1; ifile<=$nfiles; ifile++))
do
  file=$(awk <$ls_file_in "NR==$ifile {print \$1}")
  echo File nb $ifile = $file
  ln -s $file tmp_file$ifile.dat
done <$ls_file_in

# create and save 3rd-order tensor with Matlab
matlab -nosplash -nodisplay -nodesktop < make_time_lapse_input_file.m > tmp_log.out
mv time-lapse_input.mat $file_out

# clean
rm tmp*
