%
% Function FORWARD_PARAMETERS: defines input variables for the forward problem.
%
% All input.* variables should be defined here,
% to avoid hard-coding in other routines.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function [input]=forward_parameters(input)


disp(' ')
disp('------------------------------')
disp('    ENTER FORWARD_PARAMETERS  ')
disp(' defines simulation variables ')
disp('------------------------------')
disp(' ')


%------------------%
%    FILE NAMES    %
%------------------%

%% MESH
input.read_mesh=0;   % 0 -> create mesh based on acquisition geometry (create_mesh3.m)
                     % 1 -> read previously-generated mesh in input file
% file mesh_in is used if input.read_mesh==1
input.file_mesh_in='data/does-not-exist';
% file mesh_out is written down anyway
input.file_mesh_out=['mesh_sandbox_forward_4-cells-per-elec-spacing_51-layers.mat'];


%% INPUT MODEL
input.cmplx_format=1;        % I/0 format for models, 1->real/imag parts, 2->amp/phase
input.interp_model_flag=1;   % 1->interpolate input model on forward mesh, other->don't (input model has to match the mesh, then)
input.binary_model_flag=1;   % 1->round model such that it contains only 2 types of media (min/max of both real and imag. parts, specific to FL's case)
input.mean_res_flag=1;       % 1->consider background value (ie 1st cell's value) as mean value, other->take mean value of entire model
                             %    /!\ The mean value has to be accurate for accurate solutions.

% file_model_in contains the model for simulation and the initial model in case of inversion
input.file_model_in=['demodata/sandbox_forward_f1Hz/model_sandbox_true_f1Hz_251x601_h0.002.dat'];

% file_model_out contains the model that has been interpolated on the forward grid
% if model is not interpolated (input.interp_model==0), then model_out is identical to model_in
input.file_model_out=['model_sandbox_true_f1Hz_interp-on-mesh-fwd_4cpe_51l.dat'];


%% DATA
input.read_acqui=0;   % 0 -> create acquisition (ie combinations of quadripoles) based on input electrode locations
                      % 1 -> read acquisition configuration in input file
% input.mes_in is used to read acqui if read_acqui==1 
% and to read observed data in case of inversion
input.mes_in='data/does-not-exist.dat';
% file acqui_out is written if read_acqui==0
input.file_acqui_out='acqui_sandbox_24-electrodes_DP-DP.dat';

% input.file_data_out will contain the calculated synthetic data
input.file_data_out=['data-obs_sandbox_forward_f1Hz.dat'];


%--------------------%
%   MODEL GEOMETRY   %
%--------------------%

% geometry of input model
input.n1=251;        % nb of rows (nz)
input.n2=601;        % nb of columns (nx)
input.hstep=0.002;   % grid interval
  zmax=(input.n1-1)*input.hstep;

% geometry of created mesh
input.ncells_per_elec=1;   % nb of cells per electrode interval (set >1 to refine mesh laterally)
input.depth_n=1;      %     0-> generates mesh automatically with irregular depth intervals
if input.depth_n~=0   % not 0-> generate mesh with user-defined regular depth intervals
   nlayers=51
   dz_layers=zmax/(nlayers-1);
   input.depth_n=[0:nlayers-1]*dz_layers;
end

% PML
input.user_boundary_flag=1;   % 1->use extra-domain width defined below,
                              % 2->use a variation of Marios' scale,
                              % else->use Marios' scale    (see src/src_ip4di/create_mesh3.m for details)
input.extra_width=0.5;        % extra-domain width (outside domain of interest, serves as absorbing region)
input.x_margin=0;             % margin before 1st and after last electrode (still domain of interest)
                              % NB: should be a multiple of input.electrode_spacing

% max. size of element in mesh
input.hmax=[];   %input.electrode_spacing*2;


%-------------------------------%
%   ACQUISITION CONFIGURATION   %
%-------------------------------%

input.array_type=3;
% 1 -> Wenner
% 2 -> Schlumberger
% 3 -> Dipole-Dipole
% 4 -> Pole-Dipole
% 5 -> Pole-Pole

input.nb_electrodes=24;
input.electrode_spacing=0.05;
input.xelec_min=0.025;

% eventually restrict the number of quadripoles...
% (don't consider all electrode combinations, see data_2d.m for details).
input.xN_max=input.nb_electrodes;
input.xA_max=input.nb_electrodes;


% DATA TYPE
input.data_type=1;   %1->resistivities, 2->resistances, converted to resistivities in create_mesh3.m

input.sip_flag=1;
input.dc_flag=0;
input.time_lapse_flag=0;
input.ip_flag=0;
input.res2d_flag=0;


% compute Jacobian (2->yes, other->no)
input.jacobian_flag=0;


% DEBUG
input.debug=0;                    % 1->more verbose for debug
input.inv_flag=0;                 % init. inversion flag (such that we know if inversion is run or not)
input.first_time_call=1;          % init. flag (FL: first_time_call=0 is not used, anyway)
input.mesh_option_output=false;   % intern option to create_mesh3 to generate verbose output when generating mesh
input.accuracy_threshold=1e-6;    % threshold for accuracy in space locations (required accuracy depends on the dimensions of the model)

end   %end function

