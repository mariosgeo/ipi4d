To create the unique time-lapse data file required for time-lapse
inversion, first list your individual data files:
$ ls data-obs_*.dat >ls_data_files.txt

Then run
$ ./01_make_time_lapse_data_file.sh
which makes use of make_time_lapse_input_file.m to create a 3rd-order
tensor, as used in the inversion, and save it into a Mat file.

Proceed identically to generate time-lapse initial models.

Francois Lavoue
18 April 2016
