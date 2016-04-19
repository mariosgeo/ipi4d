%
% EXEC_INVERSION: MAIN PROGRAM FOR THE INVERSION OF DC, IP OR SIP DATA
%
% This program executes M. Karaoulis' IP4DI routines 
% without using the graphical user interface.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: v1, October 15, 2015.

clear all
close all


disp(' ')
disp('==============================')
disp('=    IP4DI: RUN INVERSION    =')
disp('==============================')
disp(' ')


% PATH TO SUBROUTINES
define_paths;


% DEFINE INPUT PARAMETERS
input=[];   %init. input structure
input=forward_parameters(input);
input=inversion_parameters(input);
input=plot_parameters(input);

% check inconsistencies in input parameters
abort_flag=check_inputs(input);
if(abort_flag==1); return; end   %terminate program


% READ ACQUI AND DATA
if input.time_lapse_flag==0
   input=read_data(input);
else
% time-lapse inversion
   [input]=read_4d_data(input);
end


% CREATE MESH
mesh=create_mesh3(input);


% DEFINE MEAN AND BACKGROUND VALUES
mesh=define_mean_res(input,mesh);


% DEFINE INITIAL MODEL
[input,mesh]=edit_model(input,mesh);


if input.sensitivity_analysis_flag==1
% PERFORM SENSITIVITY ANALYSIS
   [input,mesh,fem,final]=sensitivity_analysis(input,mesh);

else
% RUN INVERSION
   [fem,final]=main(input,mesh);
end

