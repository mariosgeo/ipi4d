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
dir='results_IGI_v4_gradient-TI_ACB-ok'

% choose iteration to display
itr=5;

% load path to find function plot_model_on_forward_mesh
addpath('src/src_nogui_ip4di/');

% load data
input=importdata([dir '/struct_input.dat']);
mesh=importdata([dir '/struct_mesh.dat']);
%fem=importdata([dir '/struct_fem.dat']);
final=importdata([dir '/struct_final.dat']);

% extract model at it. nb
model=final.res_param1_vs_it(:,itr);

% plot
input.plot_options.label_title=['Iteration nb ' num2str(itr)];

% plot amplitude
input.plot_options.cmplx_flag=3;
plot_model_on_forward_mesh(input,mesh,model,input.plot_options);

% plot phase
input.plot_options.cmplx_flag=4;
plot_model_on_forward_mesh(input,mesh,model,input.plot_options);

