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
%
% Update, 7 Oct. 2016
% Tends to be obsolete: rather use plot_inversion_results_from_structures.m instead.

clear all
close all

% choose results directory
dir='results/inv3_mesh-v2-15lay-807el_lagrn0.1_no-ACB_IGI-p01-plof4_force-st0/'

file_in=[dir 'interpolated_model_inv3_mesh-v2-15lay-807el_no-ACB_lagrn0.1_IGI-p01-plof4_no-force-st0.dat'];
file_in='data/training-image_profile1_Solfatara_n826x639_d0.7879x1.001896_xmin-3.018838_undersampled-x2_topo_xzv.dat'

% data format
data_format=3;   % 1 -> (x,y,m0,m1,mi,...), use itr to extract model
                 % 2 -> (x,y) (e.g. ellipses coordinates)
                 % 3 -> (x,y,z) (e.g. final model or TI)

% choose iteration to display
itr=4;

% load path to find subroutines
define_paths;

% define input variables
input=[];
input=forward_parameters(input);
input=inversion_parameters(input);
input=plot_parameters(input);

% correct for topo?
input.plot_options.plot_topo=0;

% set plot options
input.plot_options.plot_TI=1;
input.plot_options.plot_log=0;

input.plot_options.label_title=['Iteration nb ' num2str(itr)];
input.plot_options.label_title='';

input.plot_options.caxis_amp=[1 256];
input.plot_options.axis_tics_amp=[1,2,4,8,16,32,64,128,256];

input.plot_options.caxis_amp=[-0.5 0.5];
input.plot_options.axis_tics_amp=[];

input.plot_options.cmap=flipud(gray);

input.plot_options.reverse_yaxis=0;
input.plot_options.z_min=0;
input.plot_options.z_max=200;

if input.plot_options.plot_topo==1
%- mesh structure is needed to correct for topo
   mesh=importdata([dir 'mesh_struct.mat']);
end

% specify some variables that are not in current version of input and mesh yet
input.zmin_TI=0.;
mesh.dx=5.6480;
mesh.dy=5;

% load data
data=load(file_in);

% extract axis
vx=data(:,1);
vz=data(:,2);

% Solfatara: depths correspond to negative values
%vz=-vz;


if data_format==1
%- extract model at it. nb
   if input.dc_flag==1
      iitr=itr+2+1;   % convert itr to index nb (+2 because of x,z-col. before, +1 to skip it0)
      model=data(:,iitr);

   elseif input.sip_flag==1
      model_r=data(:,2*itr+1+1);   % real part of model at it. (+1 to skip it0)
      model_i=data(:,2*itr+1+2);   % imag part of model at it.
      model=model_r+1i*model_i;           % cmplx model
   end

elseif data_format==2
%- no model
   model=[];

elseif data_format==3
%- read final model
   if input.dc_flag==1
      model=data(:,3);
   elseif input.sip_flag==1
      model_r=data(:,3);   % real part
      model_i=data(:,4);   % imag part
      model=model_r+1i*model_i;   % cmplx model
   end
end


if input.dc_flag==1
   % plot resistivity
   input.plot_options.cmplx_flag=0;
elseif input.sip_flag==1
   % plot amplitude
   input.plot_options.cmplx_flag=3;
   % plot phase
   input.plot_options.cmplx_flag=4;
end


%if input.plot_options.plot_topo==1
%- mesh structure is needed to correct for topo
   plot_model_from_xzv_struct(input,mesh,vx,vz,model);

%else
%%- plot discretized or interpolated model without using Matlab structures
%   plot_model_from_xzv_nostruct(input,vx,vz,model);
%end

