%
% FUNCTION PLOT_PARAMETERS: defines plot options
%
% All input.* variables concerning plots should be
% defined here, to avoid hard-coding in other routines.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function [input]=plot_parameters(input)

%---------------------%
%    FORWARD PLOTS    %
%---------------------%

% choose plot to show (1->yes, 0->no)
input.plot_acqui=0;       % plot electrode location
input.plot_full_acqui=0;  % plot locations of all electrode combinations
input.plot_mesh=0;
input.plot_model=0;       % plot model at various stage: input, interpolation, on mesh...
input.plot_pseudosection=0;
input.plot_jacobian=0;    % plot Jacobian movie wrt quadripole combinations
input.plot_resolution=0;


%-----------------------%
%    INVERSION PLOTS    %
%-----------------------%

% choose plot to show (1->yes, 0->no)
input.log_plot=0;
input.acb_plot=0;
input.electrode_plot=0;
input.resolution_plot=0;
input.parameter_plot=0;
input.fem_plot=0;
input.pseudo_plot=0;
input.plot_ellipses=0;
input.plot_model_vs_it=0;
input.plot_interpolated_model=0;
input.plot_covariance_matrices=0;


%-------------------------------%
% OPTIONS COMMON TO FWD AND INV %
%-------------------------------%

input.plot_options.interp=0;   % 0->plot discretized model, 1->plot model after Matlab interpolation, 2->plot model after image-guided interpolation

input.plot_options.x_min=-0.025;
input.plot_options.x_max=1.225;
input.plot_options.z_min=0;
input.plot_options.z_max=0.5;

input.plot_options.caxis_amp=[275 375];   %values based on true model at 1 Hz
input.plot_options.axis_tics_amp=[275:25:375];
%input.plot_options.caxis=[345 375];   %values that enhance incomplete inversion results
%input.plot_options.ytics=[345:10:375];

input.plot_options.caxis_phi=[0.4 1.8];   %values based on true model at 1 Hz
input.plot_options.axis_tics_phi=[0.5:0.25:1.75];
%input.plot_options.caxis=[1.4 1.8];   %values that enhance incomplete inversion results
%input.plot_options.ytics=[1.4:0.1:1.8];


% color maps
map_res2d=Res2Dinv_colormap();   % Res2Dinv color map
map_GMT=cptcmap('GMT_seis');     % GMT color map

%input.plot_options.cmap=map_res2d(1:17,:);     % use Res2Dinv color map
input.plot_options.cmap=map_GMT(end:-1:1,:);   % OR use GMT map (Marios' ADD IN: reverse color map)
%input.plot_options.cmap='gray';                 % OR choose a Matlab color map ('jet', 'gray', ...)


end   %end function


