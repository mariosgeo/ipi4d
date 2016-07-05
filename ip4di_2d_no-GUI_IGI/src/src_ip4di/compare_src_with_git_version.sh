#!/bin/bash

dir1=.
dir2=../../../00_ip4di_2d_Deqiang
#$HOME/git/ip4di/ip4di_2d

ls $dir1/*.m >tmp_ls1.txt
ls $dir2/*.m >tmp_ls2.txt

# compare list of files
nf=$(wc <tmp_ls1.txt | awk '{print $1}')
echo NB OF FILES = $nf
nf=$(wc <tmp_ls2.txt | awk '{print $1}')
echo NB OF FILES = $nf

# compare files
for ((if=1; if<=$nf; if++))
do
  file1=$(awk <tmp_ls1.txt "NR==$if {print \$1}")
  file2=$(awk <tmp_ls2.txt "NR==$if {print \$1}")
  diff -q $file2 $file1
done

# clean
rm tmp*
