%
% ROUTINE INVERSION_PARAMETERS: DEFINES INPUT VARIABLES FOR THE INVERSE PROBLEM
%
% All input.* variables should be defined here,
% to avoid hard-coding in other routines.
%
% Derived from M. Karaoulis GUI routine inversion_parameters.m.
%
% Author: Francois Lavoue', Colorado School of Mines
% Last updated: March 2, 2016

function [input]=inversion_parameters(input)

disp(' ')
disp('----------------------------')
disp(' ENTER INVERSION_PARAMETERS ')
disp(' define inversion variables ')
disp('----------------------------')
disp(' ')


% Perform sensitivity analysis?
input.sensitivity_analysis_flag=0;   % 1->yes, 0->no
%if yes, user must edit file 'sensitivity_analysis.m'
%to choose the varying parameter and its range.


%---------------------%
%     FILE NAMES      %
%---------------------%

% file name for info vs. it
input.file_info='info_vs_it.txt';

% file names for final data and models
input.file_data_inv_out='data_final.dat';
input.file_model_out='model_final.dat';
input.file_model_interp='model_final_interp.dat';

% file names for intermediate data and models
input.print_intermediate_results=1;   % 1->yes, 0->no
input.file_data_inv_inter='data_inter.dat';
input.file_model_inter='model_inter.dat';

% file names for final Matlab structures,
% enabling to find any variable again
input.output_variables=1;   % 1->save structures, 0->don't to avoid storing large files
input.file_input_struct='struct_input.dat';
input.file_mesh_struct='struct_mesh.dat';
input.file_fem_struct='struct_fem.dat';
input.file_final_struct='struct_final.dat';



%------------------------%
%  TIME-LAPSE INVERSION  %
%------------------------%
input.time_lapse_flag=1;   % 1->yes, 0->no
%FL: not implemented yet in no-GUI version,
%    may cause run-time errors...

% regularize wrt time (or frequency if considering SIP data)
input.gamma=0.01;               % Lagrangian multiplier related to time regularization
input.time_smoothing_order=2;   % order of differential operator used for smoothing

% Use Active Time-Constrained regularization (ATC)
input.atc_flag=0;   % 1->yes, 0->no
%FL: not implemented yet in no-GUI version,
%    input.atc_flag=1 will cause program termination.


%-------------------------%
%   OPTIMIZATION SCHEME   %
%-------------------------%

input.inv_flag=2;
if(input.inv_flag==1); input.inv_name='Occam'; end        % 1 -> Occam (default) (suggested for a smoother model in case of noisy data)
if(input.inv_flag==2); input.inv_name='GN'; end           % 2 -> Gauss-Newton    (idem)
if(input.inv_flag==3); input.inv_name='L-M'; end          % 3 -> Levenberg-Marquardt (suggested in all other cases)
if(input.inv_flag==4); input.inv_name='Occam diff'; end   % 4 -> Occam Difference (a background model is required)
if(input.inv_flag==5); input.inv_name='GN diff'; end      % 5 -> GN Difference    (idem)
if(input.inv_flag==6); input.inv_name='??'; end           % 6 -> ??
if(input.inv_flag==7); input.inv_name='??'; end           % 6 -> ??

% GAUSS vs. QUASI NEWTON 
input.jacobian_flag=2;
% 2-> Full Jacobian calculations, provide accurate solution but time consuming.
% 1-> Quasi newton calculates Jacobian matrix once. Fast but inaccurate solution.

% REGULARISATION WEIGHT
input.decrease_lagrn_flag=0;            % decrease lagrange multiplier during inversion (1->yes, 0->no)
input.decrease_lagrn_reduction_flag=0;  % decrease lagrange multiplier reduction rate during inversion (1->yes (follow Marios' rule), 0->no)
input.lagrn_reduction=2;                % reduction rate of lambda parameter during inversion (if input.decrease_lagrn_flag==1)
input.lagrn=0.01;                       % initial value of Lagrangian mutiplier

% ACB = Active Constrained Balancing (Yi et al, 2003)
% Assign Lagrangian between two limits values based on the resolutions of the parameters.
% Use this choice in cases of borehole-data.
input.acb_flag=1;

% if acb_flag==1, specify range for Lagrangian multiplier
input.lagrn_min=0.01;
input.lagrn_max=1;    % 1 is RECOMMENDED to keep the ACB covariance matrix normalized
                      % (global weighting in ensured by input.lagrn)

% BOUNDS CONSTRAINTS: limit resistivity values
% Use this choice in cases when extremely low or high resistivity values are found.
input.limit_res=0;
input.min_res=0;
input.max_res=1.E8;


%------------------------%
%   STOPPING CRITERION   %
%------------------------%

% choice of stopping criterion
input.stop_crit=1;
% 1 -> stopping criterion considers the total objective function that is effectively minimized (RECOMMENDED)
%          ie C(m) = C_D(m) + lagrn*C_M(m)
%      with C_D(m) = (log(dobs)-log(dcal(m)))' Wd (log(dobs)-log(dcal(m)))
%       and C_M(m) = log(m)' Cm^-1 log(m)
% 2 -> stopping criterion considers only the data term C_D(m) of the minimized objective function
% 3 -> stopping criterion considers a normalized data misfit, computed as
%         nRMS = sum_i=1:ndata |dobs_i-dcal_i(m)|^2 / |dobs_i|^2
%              = (dobs-dcal(m))' Wd (dobs-dcal(m))
%      with Wd = diag( 1/|dobs_i|^2 )
% 4 -> stopping criterion considers a weighted data misfit, computed as
%         wRMS = (dobs-dcal(m))' Wd (dobs-dcal(m))
%      with user-defined Wd (=Identity in current version of the code)

% minimum convergence rate
% (in %: inversion stops if the relative decrease of the objective
%  function between two iterations is smaller than this rate)
input.conv_rate=1;

% max nb of iterations
input.itn=100;


%--------------------%
%  OTHER PARAMETERS  %
%--------------------%

% background model
input.bgr_res_flag=0;
if input.bgr_res_flag==1 || input.bgr_res_flag==2
   input.par_nm='background_file';
end

% boundary conditions
mesh.bc_option=2;   % 1 -> Dirichlet, 2 -> mixed boundary conditions

% ?? parameters
input.current=1;       %?? not used...
input.strike_dist=0;   %?? not used...



%--------------------------%
%  IMAGE-GUIDED INVERSION  %
%--------------------------%

input.image_guidance=2;
% 0-> don't use image-guidance
% 1-> use initcm1.m and 4-dimensional smoothing (Zhou et al., 2014)
% 2-> use initcm2.m and directional Laplacian filters (Lavoue et al., in prep)
% 3-> read covariance matrix in a file

% image-guided interpolation of final model (1->yes, 0->no)
input.image_guided_interpolation=0;

% if input.image_guidance==1 or 2, use training image
% size of training image
% /!\ SU conventions
input.n1_TI=251;   %nz
input.n2_TI=601;   %nx
input.hstep=0.002;

% file name for training image
input.training_image=['demodata/sandbox_inversion_f1Hz/image_sandbox_251x601_h0.002.bin'];

% if input.image_guidance==3, read covariance matrix in a file
input.file_covariance='data/does-not-exist';

% IMAGE GUIDANCE PARAMETERS
mesh.scale=1;   %FL: should probably not be used anymore (ie let to 1),
                %    because it is redundant with the spatially varying Lagrangian distribution defined by ACB.

input.IGI.p_lof1=10;   % half-width of isotropic Gaussian window for LocalOrientFilter
input.IGI.p_lof2=10;   % (in nb of cells, should be at least 4 according to JZ, 
                       %  but may also depend on the nb of TI cells in 1 inversion cell)
                       % NB:it is also the diffusivity alpha used for interpolation (see DH's CWP report)

%NB: the half-widths p_lofs relate to a second smoothing applied after building the tensor T (K_rho in Weckert, )
%FL: according to Elias Arias' experience on dip filtering, an anisotropic Gaussian
%    window is suitable if the training image contains anisotropic structures.
%    If the structures display more variations in direction 1 (vertical), we may then
%    use a larger p_lof1 than p_lof2 (e.g. [p_lof1,p_lof1/2] or [p_lof1,p_lof1/3]).
%      Example:
%      anisotropic_ratio=3;
%      input.IGI.p_lof1=round(input.n1_TI/mesh.m2)*anisotropic_ratio;
%      input.IGI.p_lof2=round(input.n2_TI/mesh.m1);

input.IGI.p0=0.0;    %JZ: These are the p0, p1 in my paper.
input.IGI.p1=1.0;    %    The larger the second parameter, the more anisotropic the tensors,
                     %    can be set to 3 or 4, but for geological cross-section image please leave as 1.0.
%DH's doc for invertStructure:
% p0 emphasizes overall amplitude and p1 emphasizes linearity (ie anisotropy).
% For amplitude-independent tensors with all eigenvalues av equal to one, set p0 = 0.0.
% To enhance linearity, set p1 > 1.0. To simply invert (and normalize) these tensors, set p0 = p1 = 1.0.

%% Parameters of local smoothing filters (lsf) for local semblance computation
input.IGI.p_lsf1_1=input.IGI.p_lof1;   % lsf1 is related to coherence (not used)
input.IGI.p_lsf1_2=input.IGI.p_lof2;

input.IGI.p_lsf2_1=input.IGI.p_lof1;   % half-width of 1st smoothing filter for lsf2 (related to semblance, used)
input.IGI.p_lsf2_2=input.IGI.p_lof2;   % half-width of 2nd smoothing filter


end   %end function


