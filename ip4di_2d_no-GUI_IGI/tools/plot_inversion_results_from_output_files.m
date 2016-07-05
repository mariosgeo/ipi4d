%
% MATLAB PROGRAM TO PLOT INVERSION RESULTS
% from output ASCII files.
%
% NB:
% - works best with Java 1.6 (while IGI needs Java 1.7).
% - calls the input functions 'forward_parameters.m',
%   'inversion_parameters.m' and 'plot_parameters.m'
%   which may be adapted accordingly (in particular, make
%   sure they are consistent with the settings used for the
%   inversion that have generated the model to be plotted).
%
% F. Lavoue', Colorado School of Mines
% October 25, 2015

clear all
%close all

%= USER PARAMETERS =%
% choose results directory
dir='results_IGI_v5_true-TI_ACB-ok_DIVERGENCE'

% choose iteration to display
itr=1;

% define input variables
input=[];
input=forward_parameters(input);
input=inversion_parameters(input);
input=plot_parameters(input);

% plot options
input.plot_options.label_title=['Iteration nb ' num2str(itr)];
%= END USER PARAMETERS =%


% load path to find function plot_model
addpath('src/src_nogui_ip4di/');
addpath('src/src_ip4di/');
addpath('cmaps')

% load data
model_vs_it=load([dir '/rhoc_ri_sandbox_inv_inter.mod']);
%FL DEBUG
model_vs_it=load('data-obs_sandbox_mesh-inv/rhoc_ri_sandbox_mesh-inv_true_f1Hz.mod');

% extract model at it. nb
vx=model_vs_it(:,1);
vz=model_vs_it(:,2);

if input.dc_flag==1
   model=model_vs_it(:,itr+2);   % +2 because of x,z-col. before

elseif input.sip_flag==1
   model_r=model_vs_it(:,2*itr+1);   % real part of model at it.
   model_i=model_vs_it(:,2*itr+2);   % imag part of model at it.
   model=model_r+1i*model_i;         % cmplx model
end

% plot amplitude
input.plot_options.cmplx_flag=3;
plot_model_from_xzv(input,vx,vz,model);

% plot phase
if input.sip_flag==1
   input.plot_options.cmplx_flag=4;
   plot_model_from_xzv(input,vx,vz,model);
end

