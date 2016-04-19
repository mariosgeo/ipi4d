% 	  /* ******************************************************************
% 	   *                             RMS( )                               *
% 	   *------------------------------------------------------------------*
% 	   * This function finds the RMS error and ends the procedure         *
% 	   ********************************************************************/
%
% Author: Marios Karaoulis, Colorado School of Mines
%
% Modified by Francois Lavoue, April 2016
%   Modifications include
%   1) some explicit comments for more readability.
%   2) the definition of various measures of misfits:
%      a) the normalized RMS already defined by Marios 
%         (but it is not the objective function that we actually minimize in invert_cntr.m)
%      b) the total objective function (actually minimized in invert_cntr.m)
%             C(m) = C_D(m) + lambda*C_M(m)
%         with the data term C_D(m) = || log(dobs) - log(dcal(m)) ||^2
%         and the model term C_M(m) = m' Cm^-1 m
%   3) some subsequent change of notations:
%      - 'tsum' is now 'nrms' (normalized RMS),
%      - 'fem.rms_sum1' is now 'fem.nrms',
%      - 'fem.rms_sum2' is now 'fem.nrms2' (normalized RMS at previous iteration, not used anymore, use fem.rms_crit2 instead),
%   4) the definition of the locally-weighted ACB covariance matrix here rather than in
%      invert_cntr.m because we need it to compute the model term of the objective function.

function [ex,fem,mesh]=rms_4d(itr,ip_cnt,input,mesh,fem)

disp(' ')
disp('---------------------------------------')
disp('             ENTER RMS_4D              ')
disp(' computes RMS for time-lapse inversion ')
disp('---------------------------------------')
disp(' ')

% init. misfits for current iteration
nrms=0;
fem.rms=0;         % simple RMS misfit
fem.nrms=0;        % normalized RMS
fem.wrms=0;        % weighted RMS
fem.mrms=0;        % model misfit
fem.obfj=0;        % total objective function
fem.obfj_data=0;   % data term of objective function
fem.obfj_model=0;  % model term of objective function
fem.rms_crit=0;    % criterion considered for stopping optimization
                   %   (can be either RMS, nRMS, wRMS, data or total objective function,
                   %    depending on user-defined input.stop_crit)

% init. exit signal
% (0->CONTINUE, 1->STOP, 2->invert last IP component and STOP)
ex=0;


% recast model and data in single vectors
tmp_d4_obs=[];
tmp_d4_cal=[];
tmp_d4_mod=[];
for i=1:input.num_files;
    tmp_d4_obs = [tmp_d4_obs ; input.d4_real_data(:,i)];
    tmp_d4_cal = [tmp_d4_cal ; fem.d4_array_model_data(:,i)];
    tmp_d4_mod = [tmp_d4_mod ; mesh.d4_res_param1(:,i)];
end
        

%------------------------------------%
%        COMPUTE DATA MISFIT         %
%------------------------------------%

% Find RMS error
for i=1:input.num_mes*input.num_files
    %compute normalized RMS of data misfit
    % C = sum_i=1:Nd  | dobs_i - dcal_i(m) |^2 / | dobs_i |^2
    %FL: This is NOT the objective function which is minimized.
    %    In this formulation, the residuals are weighted by the inverse of the norm of the observed data.
    %    As a consequence, all residuals have a similar impact on the misfit value (since low amplitudes are re-inflated).
    nrms = nrms + ( (tmp_d4_obs(i)-tmp_d4_cal(i))' * (tmp_d4_obs(i)-tmp_d4_cal(i)) ) ./ (tmp_d4_obs(i)'*tmp_d4_obs(i));
end

% compute RMS of data misfit (without weight/normalization)
% C = || dobs - dcal(m) ||^2
fem.rms = (tmp_d4_obs-tmp_d4_cal)' * (tmp_d4_obs-tmp_d4_cal);
% convert RMS in percentage of error
fem.rms = sqrt(fem.rms/(input.num_mes*input.num_files))*100;
disp(sprintf('RMS = %f %%',fem.rms));

% convert normalized RMS in percentage of error
fem.nrms = sqrt(nrms/(input.num_mes*input.num_files))*100;
disp(sprintf('Normalized RMS = %f %%',fem.nrms));

%compute weighted RMS of data misfit (with data weighting matrix Wd)
% C = (dobs-dcal(m))' Wd (dobs-dcal(m))
fem.wrms = (tmp_d4_obs-tmp_d4_cal)' * input.Wd_d4 * (tmp_d4_obs-tmp_d4_cal);
% convert weighted RMS in percentage of error
fem.wrms = sqrt(fem.wrms/(input.num_mes*input.num_files))*100;
disp(sprintf('Weighted RMS = %f %%',fem.wrms));

%compute RMS of log misfit, as used in invert_cntr.m for inversion
% C = [ log(dobs)-log(dcal(m)) ]' Wd [ log(dobs)-log(dcal(m)) ]
fem.objf_data = fem.e' * input.Wd_d4 * fem.e;

% display log RMS and percentage of error
% (does not convert fem.objf_data itself to preserve the quantity that is effectively minimized)
disp(sprintf('LOG RMS = %f',fem.objf_data));
disp(sprintf('LOG RMS = %f %%',sqrt(fem.objf_data/(input.num_mes*input.num_files))*100));



%------------------------------------%
%        COMPUTE MODEL TERM          %
%------------------------------------%

% C_M(m) = m' Cm^-1 m
fem.mrms = tmp_d4_mod' * fem.CC * tmp_d4_mod;
disp(sprintf('C_M(m) = %f',fem.mrms));
disp(sprintf('C_M(m) = %f %%',sqrt(fem.mrms/mesh.num_param)*100));

% C_M(m) = log(m)' Cm^-1 log(m)
fem.objf_model = log10(tmp_d4_mod)' * fem.CC * log10(tmp_d4_mod);
disp(sprintf('C_M(log m) = %f',fem.objf_model));
disp(sprintf('C_M(log m) = %f %%',sqrt(fem.objf_model/mesh.num_param)*100));


%--------------------------------------%
%   COMPUTE TOTAL OBJECTIVE FUNCTION   %
%--------------------------------------%

% C(m) = C_D(m) + lambda*C_M(m)
fem.objf = fem.objf_data + input.lagrn*fem.objf_model;
disp(' ')
disp('---------------')
disp(sprintf('C_D(m) = %f',fem.objf_data));
disp(sprintf('C_M(m) = %f',fem.objf_model));
disp(sprintf('lambda*C_M(m) = %f',input.lagrn*fem.objf_model));
disp(sprintf('lambda*C_M(m) / C_D(m) = %f',input.lagrn*fem.objf_model/fem.objf_data));
disp(sprintf('C(m) = %f',fem.objf));
disp('---------------')
disp(' ')


%
% CHOOSE STOPPING CRITERION
%
if input.stop_crit==1
   %use total objective function
   % (RECOMMENDED because it is the one we minimize in invert_cntr.m)
   fem.rms_crit=fem.objf;
elseif input.stop_crit==2
   %use only data term of objective function
   fem.rms_crit=fem.objf_data;
elseif input.stop_crit==3
   %use normalized RMS
   fem.rms_crit=fem.nrms;
elseif input.stop_crit==4
   %use weighted RMS
   fem.rms_crit=fem.wrms;
end


% Update for first iteration
if (itr==1)
   fem.rms_crit2=fem.rms_crit;
end


% Check divergence
% (should not occur if the total objective function has been chosen as criterion, input.stop_crit=1,
%  otherwise, it means there is a bug in the model update in invert_cntr.m.)
if (fem.rms_crit>fem.rms_crit2 && itr>1)
    disp(sprintf('Divergence occurred: crit(it-1) = %f, crit(it) = %f',fem.rms_crit2,fem.rms_crit));
    disp(sprintf('                     nRMS(it-1) = %f, nRMS(it) = %f',fem.nrms,fem.nrms));
    fem.rms_crit=fem.rms_crit2;   %keep previous value as final one

    if input.ip_flag==1 && ip_cnt==1
    % If we have IP data, update everything to continue with IP inversion 
        ex=2;

    elseif input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==2)
    %else, exit signal 1 ==> program termination
        ex=1;
    end
end


% Compute relative misfit decrease since last iteration
rel_dec=abs((fem.rms_crit2-fem.rms_crit)/fem.rms_crit2);

% Check convergence rate
if ( rel_dec<(input.conv_rate/100) && itr>1)
    disp(sprintf('Convergence is slower than user-defined rate of %f%%: crit(it-1) = %f, crit(it) = %f',...
                  input.conv_rate,fem.rms_crit2,fem.rms_crit));
    disp(sprintf('                                                      nRMS(it-1) = %f, nRMS(it) = %f',...
                  fem.nrms,fem.nrms));
    fem.rms_crit=fem.rms_crit2;

    if input.ip_flag==1 && ip_cnt==1
    % If we have IP data, update everything to continue with IP inversion 
        ex=2;
    elseif input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==2)
    %else, exit signal 1 ==> program termination
        ex=1;
    end
end


if ex==0
    %before going to next iteration,
    %update RMS, data and model of 'previous' iteration
    %(don't do it for ex==2, in order to invert the 2nd IP component
    % with the current model and data)
    mesh.d4_res_param2=mesh.d4_res_param1;
    fem.nrms2=fem.nrms;
    fem.d4_array_model_data2=fem.d4_array_model_data;
end


if itr==input.itn+1 && ex==0
%if max. nb of iterations has been reached, terminate program
%(unless if there is still an IP component to invert)

    if input.ip_flag==1 && ip_cnt==1
    % If we have IP data update everything to continue with IP inversion
       ex=2;
       disp('********IP INVERSION STARTS**********');

    elseif input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==2)
    % otherwise, terminate program
       ex=1;
       disp('********PROGRAM TERMINATION**********');
    end
end


end   %end function rms_4d

