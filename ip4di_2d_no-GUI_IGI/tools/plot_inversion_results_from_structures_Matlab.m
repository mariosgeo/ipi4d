%
% MATLAB PROGRAM TO PLOT INVERSION RESULTS
% from saved Matlab structures.
%
% NB: - works best with Java 1.6
%       (while IGI needs Java 1.7)
%
% F. Lavoue', Colorado School of Mines
% October 24, 2015

clear all
close all

% choose results directory
dir='results/test9_inversion_v3-Lcurve-807R-xas/'

% choose iteration to display
itr=4

% load path to find function plot_model_on_forward_mesh
addpath('src/src_nogui_ip4di/');

% load data
input=importdata([dir 'input_struct_lagrn0.1.mat']);
mesh=importdata([dir 'mesh_struct_lagrn0.1.mat']);
%fem=importdata([dir '/fem_struct.mat']);
final=importdata([dir 'final_struct_lagrn0.1.mat']);

% extract model at it. nb
model=final.res_param1_vs_it(:,itr);

% plot L-curve
input.plot_options.label_title=['Iteration nb ' num2str(itr)];

% plot amplitude
input.plot_options.cmplx_flag=3;
plot_model_on_forward_mesh(input,mesh,model)

% plot phase
if input.sip_flag==1
   input.plot_options.cmplx_flag=4;
   plot_model_on_forward_mesh(input,mesh,model,input.plot_options);
end

