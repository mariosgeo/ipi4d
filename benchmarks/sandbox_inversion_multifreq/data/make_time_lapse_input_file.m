clear all
close all

% read nb of input files
nfiles=load('tmp_nfiles.txt')

% create 3rd-order tensor
for ifile=1:nfiles
    data=load(['tmp_file' num2str(ifile) '.dat']);
    time_lapse_data{ifile}=data;
end

% save 3rd-order tensor as a Mat file
save('time-lapse_input.mat','time_lapse_data');
%(needs to be read using importdata)
