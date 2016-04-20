%
% MAIN INVERSION ROUTINE
%
% Author: Marios Karaoulis, Colorado School of Mines

% Modified by Francois Lavoue', Colorado School of Mines, September 2015.
% Modifications include:
% 1) some explicit comments for more readability,
% 2) some change of notations (TO DO),
%    in particular the structure 'final' should be denoted 'inv' because it gathers all inversion variables
%              and the structure 'fem' may be denoted 'fwd' since it concerns forward variables...
% 3) some ordering changes:
%    - Lagrangian multiplier is now updated after parameter update for consistency.

function [fem,final]=main(input,mesh)

%global fem final

disp(' ')
disp('===================================')
disp('=   ENTER MAIN INVERSION PROGRAM  =')
disp('===================================')
disp(' ')

%if matlabpool('size') == 0
%% parallelize
%  matlabpool open
%  matlabpool_size=matlabpool('size')
%end

% init. forward and inverse structures        
fem=[];
final=[];

% DEFINE SMOOTHING MATRIX 
% (via classical covariance or image-guidance)
disp('******** SMOOTH MATRIX *******');
mesh=smooth_mtx_surface4(input,mesh);

% Calculate k and g
[mesh.k,mesh.g]=rwaven(mesh.probe_spacing,mesh.probe_spacing*mesh.max_n);
% rdis=double([mesh.probe_spacing:mesh.probe_spacing:mesh.max_n*mesh.probe_spacing]);
% k_init= logspace(-2,0.3,5);
% [mesh.k,mesh.g,obj,err]=xu_inversion(rdis,k_init);

disp('****** INVERSION  STARTS *****');

if input.time_lapse_flag==0
% NO TIME-LAPSE INVERSION 

    % loop over IP components
    % (ip_num=2 only in case of IP data, not SIP)
    for ip_cnt=1:input.ip_num

        % LOOP OVER ITERATIONS
        for itr=1:input.itn+1

            %start chrono
            t1=cputime;

            %compute synthetic data and Jacobian matrix
            fem=mes_control_fast(itr,input,mesh,fem,0);

            %compute misfit
            [ex,fem,mesh]=rms(itr,ip_cnt,input,mesh,fem);

            %update variables and print intermediate results
            %for current iteration (before model update)
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);

            %check exit signal (and eventually skip parameter update below)
            if(ex==1); disp('********PROGRAM TERMINATION**********'); end
            if(ex==1 || ex==2); final.itr=itr-1; break; end

            %update model parameters
            [mesh,fem,input]=invert_cntr(itr,ip_cnt,input,mesh,fem);

            %update forward model on finie-element mesh
            [fem,mesh,input]=prop_change(itr,input,mesh,fem);

            %update Lagrangian multiplier
            [input]=update_lagran(itr,ip_cnt,input);

            %plot updated model
            if(input.plot_model_vs_it==1); auto_contour(1,itr,ip_cnt,input,mesh,fem); end

            %save resolution and Jacobian matrices
            %(after model update because resolution is computed in invert_cntr.m)
            final=info_save(1,itr,ip_cnt,input,mesh,fem,final);

            %print CPU time
            t2=cputime;
            disp(sprintf('TIME FOR ITERATION = %f',t2-t1));

        end   %end loop over iterations

        if itr==input.itn+1
           final.itr=itr
        end

        if input.ip_flag==1
            [input,mesh]=ip_calc(input,mesh);
        else
            break;
        end
    end   %end loop over IP components



elseif input.time_lapse_flag==1
% TIME-LAPSE INVERSION

    % loop over IP components
    % (ip_num=2 only in case of IP data, not SIP)
    for ip_cnt=1:input.ip_num   

        % LOOP OVER ITERATIONS
        for itr=1:input.itn+1

            %compute synthetic data for each time step
            [fem,mesh]=d4_prepare_data(itr,input,mesh,fem);

            %compute misfit
            [ex,fem,mesh]=rms_4d(itr,ip_cnt,input,mesh,fem);
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);

            %check exit signal
            %(0->CONTINUE, 1->STOP, 2->invert 2nd IP component and STOP)
            if (ex==1 || ex==2)
               final.itr=itr-1;
               if(ex==1); disp('********PROGRAM TERMINATION**********'); end
               break;
            end

            %update Lagragian multiplier
            [input]=update_lagran(itr,ip_cnt,input);

            %update model parameters
            [mesh]=kim_inversion2(input,mesh,fem);

            %plot updated model
            if(input.plot_model_vs_it==1); auto_contour_4d(itr,input,mesh,fem); end
        end   %end for itr

        if itr==input.itn+1
           final.itr=itr;
        end

        if input.ip_flag==1
            [input,mesh]=ip_calc(input,mesh);
        else
            break;
        end   
    end   %end for ip_cnt
end   %end if time-lapse


% Here save outputs...
info_save(3,itr,ip_cnt,input,mesh,fem,final);
if input.output_variables==1
% save structures of variables
   save(input.file_input_struct,'input');
   save(input.file_mesh_struct,'mesh'); 
   save(input.file_final_struct,'final');
   %save(input.file_fem_struct,'fem');   %except fem because to heavy
   % (interesting variables have been saved in final anyway)
end


end   %end function main


