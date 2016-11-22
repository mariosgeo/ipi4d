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
%dir='results/Lcurve_mesh-v3-31lay-807el-dz5m_m0-rho-app10_no-std_no-ACB_IGI-plof4-p0-1_force-st0/'
%base_file_out='fig-inv1_mesh-v3-31lay-807el-dz5m_m0-rho-app10_no-std_no-ACB_IGI-plof4-p0-1_force-st0'

dir='results/Lcurve_mesh-v3-31lay-807el-dz5m_m0-rho-app10_no-std_no-ACB_no-IGI/'
base_file_out='fig-inv1_mesh-v3-31lay-807el-dz5m_m0-rho-app10_no-std_no-ACB_no-IGI'

% plot TI on top of model (does not work...)
plot_TI=0;
file_TI='data/training-image_profile1_Solfatara_n826x639_d0.7879x1.001896_xmin-3.018838_undersampled-x2_topo_xzv.dat'

% choose iteration to display
itr=4

% choose lagrn multiplier
lagrn=0.01
suffix=['_lagrn' num2str(lagrn)];
%suffix='_inter';

% load path to find subroutines
define_paths;

% choose to use Matlab or Octave
matlab_flag=1;   % 1 -> Matlab, else -> Octave (does not work yet)

%- read structure files
if matlab_flag==1
   input=importdata([dir 'input_struct' suffix '.mat']);
   mesh=importdata([dir 'mesh_struct' suffix '.mat']);
   %fem=importdata([dir '/fem_struct.mat']);
   final=importdata([dir 'final_struct' suffix '.mat']);
else   % Octave
   load([dir 'input_struct' suffix '.mat']);
end

% title
input.plot_options.label_title=['IGI (p_0=0, p_1=1, p_{lof}=4), no ACB, \lambda = ' num2str(lagrn) ', it. ' num2str(itr)];

% plot topo?
input.plot_options.plot_topo=1;

% choose type of plot
input.plot_options.interp= 2;  % -1 -> plot discretized model on forward mesh
                               %  0 -> plot discretized model without mesh
                               %  1 -> plot Matlab-interpolated model
                               %  2 -> plot IGI-interpolated model

% new sampling for Matlab interpolation
dx_int=1;
dz_int=1;

% specify some variables that are not in current version of input and mesh yet
input.zmin_TI=0.;
mesh.dx=5.6480;
mesh.dy=5;

% extract model at it. nb
iitr=itr+1;   % convert itr to index nb
model=final.res_param1_vs_it(:,iitr);

% set plot options
input.plot_options.plot_TI=0;
input.plot_options.plot_log=1;

input.plot_options.caxis_amp=[1 256];
input.plot_options.axis_tics_amp=[1,2,4,8,16,32,64,128,256];

input.plot_options.reverse_yaxis=0;
input.plot_options.z_min=-58.470;
input.plot_options.z_max=200;

input.plot_options.x_min=-10.648;
input.plot_options.x_max=643.20;

if input.dc_flag==1
   % plot resistivity
   input.plot_options.cmplx_flag=0;
elseif input.sip_flag==1
   % plot amplitude
   input.plot_options.cmplx_flag=3;
   % plot phase
   input.plot_options.cmplx_flag=4;
end


if input.plot_options.interp==-1
%- plot discretized model on forward mesh

   input.plot_options.reverse_yaxis=1;
   input.plot_options.z_min=-200;
   input.plot_options.z_max=58.470;
   suffix_out=[];

   plot_model_on_forward_mesh(input,mesh,model);

elseif input.plot_options.interp==0
%- plot discretized or interpolated model without mesh

   plot_model_from_xzv_struct(input,mesh,mesh.param_x,-mesh.param_y,model);
   suffix_out='_no-mesh';

elseif input.plot_options.interp==1
%- plot Matlab-interpolated model

   % title
   input.plot_options.label_title=[input.plot_options.label_title ', Matlab-interpolated model'];
   suffix_out='_iso-interp';

   % new grid
   vx_int=[min(mesh.param_x):dx_int:max(mesh.param_x)];
   vz_int=[min(mesh.param_y):dz_int:max(mesh.param_y)];

   % interpolate
   [MX,MZ]=meshgrid(vx_int,vz_int);
   model_interp=griddata(mesh.param_x,mesh.param_y,model,MX,MZ,'linear');

   % plot
   plot_model_from_xzv_struct(input,mesh,MX(:),MZ(:),model_interp(:));

elseif input.plot_options.interp==2
%- plot IGI-interpolated model

   % title
   input.plot_options.label_title=[input.plot_options.label_title ', IGI-interpolated model'];
   suffix_out='_IGI-interp';

   % interpolate
   [model_interp,mesh]=image_guided_interpolation(input,mesh,model);
   [MX,MZ]=meshgrid(mesh.vx_int,mesh.vz_int);

   % plot
   plot_model_from_xzv_struct(input,mesh,MX(:),MZ(:),model_interp(:));
end


%- plot TI on top of model
if plot_TI==1
   data=load(file_TI);
   vx=data(:,1);
   vz=data(:,2);
   model=data(:,3);

   input.plot_options.plot_TI=1;
   input.plot_options.plot_log=0;
   input.plot_options.caxis_amp=[-0.5 0.5];
   input.plot_options.cmap=flipud(gray);
   plot_model_from_xzv_struct(input,mesh,vx,vz,model);
end


%- suggestion for printing to EPS file
disp(' ');
disp('Resize and print to EPS file:');
file_out=['figs/' base_file_out '_lagrn' num2str(lagrn) '_it' num2str(itr) suffix_out '.eps'];
disp(sprintf('print(''-depsc'',''%s'')',file_out));
disp(' ');

