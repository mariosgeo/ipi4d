%
% FUNCTION CHECK_INPUTS: check input variables and correct them
%                        or terminate program in case of inconsistency.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function [abort_flag]=check_inputs(input)

 abort_flag=0;   %init.

 % need to compute Jacobian to plot resolution
 if input.plot_resolution==1; input.jacobian_flag=2; end

 if input.inv_flag>0 && input.atc_flag==1
    disp('ATC NOT IMPLEMENTED IN NO-GUI VERSION.');
    disp('Please set input.atc_flag=0.');
    abort_flag=1;   %terminate program
 end

 if input.inv_flag>0 && input.image_guidance==1
    disp('IGI VIA 4-DIRECTIONAL SMOOTHING NOT VALIDATED YET.');
    disp('Please set input.image_guidance=0, 2 or 3.');
    abort_flag=1;   %terminate program
 end

 if input.dc_flag==1
 %- in the DC case, just plot resistivity
    input.plot_options.cmplx_flag=0;
 end

end   %end function check_inputs

