To create the unique time-lapse data file required for time-lapse
inversion, first list your individual data files in a file:
$ ls data-obs_*.dat >ls_data_files.txt

Then run
$ ./make_time_lapse_data_file.sh
which makes use of make_time_lapse_data_file.m to create a 3rd-order
tensor, as used in the inversion, and save it into a Mat file.

Francois Lavoue
18 April 2016
