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
dir='results/test11_inversion_v5-Lcurve-807rho-app-xz-ACB/'

% choose iteration to display
itr=2

% choose lagrn multiplier
lagrn=10

% load path to find function plot_model_on_forward_mesh
addpath('src/src_nogui_ip4di/');

% load data
input=importdata([dir 'input_struct_lagrn' num2str(lagrn) '.mat']);
mesh=importdata([dir 'mesh_struct_lagrn' num2str(lagrn) '.mat']);
%fem=importdata([dir '/fem_struct.mat']);
final=importdata([dir 'final_struct_lagrn' num2str(lagrn) '.mat']);

% extract model at it. nb
model=final.res_param1_vs_it(:,itr);

% set plot options
%input.plot_options.plot_log=1;
input.plot_options.label_title=['Iteration nb ' num2str(itr)];
%input.plot_options.caxis_amp=[1 256];   %values based on true model at 1 Hz
%input.plot_options.axis_tics_amp=[1,2,4,8,16,32,64,128,256];
%input.plot_options.ymin=-200;
%input.plot_options.ymax=0;

if input.dc_flag==1
   % plot resistivity
   input.plot_options.cmplx_flag=0; plot_model_on_forward_mesh(input,mesh,model)
elseif input.sip_flag==1
   % plot amplitude
   input.plot_options.cmplx_flag=3; plot_model_on_forward_mesh(input,mesh,model)
   % plot phase
   input.plot_options.cmplx_flag=4; plot_model_on_forward_mesh(input,mesh,model)
end

