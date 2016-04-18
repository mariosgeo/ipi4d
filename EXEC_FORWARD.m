%
% EXEC_FORWARD: MAIN PROGRAM FOR THE RESOLUTION OF THE FORWARD ELECTRICAL PROBLEM
%
% This program, derived from original IP4DI forward_modelling.m,
% executes M. Karaoulis' IP4DI routines without using the graphical user interface.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: v1, October 15, 2015.

clear all
close all

global input mesh fem

disp(' ')
disp('=======================================')
disp('=    IP4DI: RUN FORWARD SIMULATION    =')
disp('=======================================')
disp(' ')


% PATH TO SUBROUTINES
define_paths;


% DEFINE INPUT PARAMETERS
input=[];   %init. input structure
input=forward_parameters(input);
%input=inversion_parameters(input);
input=plot_parameters(input);

% check inconsistencies in input parameters
abort_flag=check_inputs(input);
if abort_flag==1
   disp('CHECK INPUT PARAMETERS')
   return
end   %terminate program


% CREATE ACQUI
if input.read_acqui==0
   input=data_2d(input);
end


% READ ACQUI
input=read_data(input);


% CREATE MESH
mesh=create_mesh3(input);


% DEFINE MEAN AND BACKGROUND VALUES
mesh=define_mean_res(input,mesh);


%% DEFINE MODEL
[input,mesh]=edit_model(input,mesh);


% RUN FORWARD SOLVER
fem=forward_solver(input,mesh);


% PLOT DATA
if input.plot_pseudosection==1
   data=load(input.file_data_out);
   plot_pseudosection(input,data);
end


% PLOT JACOBIAN
if input.plot_jacobian==1
   plot_Jacobian(input,mesh,fem);
end


% RESOLUTION ANALYSIS
if input.plot_resolution==1
   plot_resolution(input,mesh,fem)
end

