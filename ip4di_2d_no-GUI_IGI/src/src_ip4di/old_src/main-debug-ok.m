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

if matlabpool('size') == 0
% parallelize
%  matlabpool open
%  matlabpool_size=matlabpool('size')
end

% init. forward and inverse structures        
fem=[];
final=[];

%FL debug
if input.debug==1
   input.debug_id=fopen('debug.out','wt');
   fprintf(input.debug_id,'   DEBUG FILE   \n\n');
   fclose(input.debug_id);
   input.debug_id=fopen('debug.out','a');   %re-open in 'append' mode
end

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

        %FL debug
        if(input.debug==1); fprintf(input.debug_id,'Enter iteration loop (nit_max=%d) \n',input.itn+1); end

        % LOOP OVER ITERATIONS
        for itr=1:input.itn+1

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'\n===============\nIteration nb=%d \n',itr); end

            %start chrono
            t1=cputime;

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Solve forward problem \n'); end

            %compute synthetic data and Jacobian matrix
            fem=mes_control_fast(itr,input,mesh,fem,0);

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Enter rms.m: compute misfit \n'); end
            if(input.debug==1); fprintf(input.debug_id,'             input.lagrn=%f \n',input.lagrn); end

            %compute misfit
            [ex,fem,mesh]=rms(itr,ip_cnt,input,mesh,fem);

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Out of rms.m: objf=%f \n',fem.objf); end
            if(input.debug==1); fprintf(input.debug_id,'              input.lagrn=%f \n',input.lagrn); end
            if(input.debug==1); fprintf(input.debug_id,'              exit=%d \n',ex); end
            if(input.debug==1); fprintf(input.debug_id,'Enter info_save.m \n'); end

            %update variables and print intermediate results
            %for current iteration (before model update)
            final=info_save(0,itr,ip_cnt,input,mesh,fem,final);

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Out of info_save.m: objf=%f \n',fem.objf); end
            if(input.debug==1); fprintf(input.debug_id,'                    input.lagrn=%f \n',input.lagrn); end

            %check exit signal (and eventually skip parameter update below)
            if(ex==1); disp('********PROGRAM TERMINATION**********'); end
            if(ex==1 || ex==2); final.itr=itr-1; break; end

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Enter invert_cntr.m and prop_change.m \n'); end

            %update model parameters
            [mesh,fem,input]=invert_cntr(itr,ip_cnt,input,mesh,fem);

            %update forward model on finie-element mesh
            [fem,mesh,input]=prop_change(itr,input,mesh,fem);

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Out of prop_change.m: objf=%f \n',fem.objf); end
            if(input.debug==1); fprintf(input.debug_id,'                      input.lagrn=%f \n',input.lagrn); end
            if(input.debug==1); fprintf(input.debug_id,'Enter update_lagran.m \n'); end

            %update Lagrangian multiplier
            [input]=update_lagran(itr,ip_cnt,input);

            %FL debug
            if(input.debug==1); fprintf(input.debug_id,'Out of update_lagran.m: objf=%f \n',fem.objf); end
            if(input.debug==1); fprintf(input.debug_id,'                        input.lagrn=%f \n',input.lagrn); end

            %plot updated model
            if input.plot_model_vs_it==1
               auto_contour(1,itr,ip_cnt,input,mesh,fem);
            end

            %print CPU time
            t2=cputime;
            disp(sprintf('TIME FOR ITERATION = %f',t2-t1));
            if(input.debug==1); fprintf(input.debug_id,'TIME FOR ITERATION = %f \n',t2-t1); end

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
            %(0->CONTINUE, 1->STOP, 2->invert last IP component and STOP)
            if (ex==1 || ex==2)
               final.itr=itr-1;
               if(ex==1); disp('********PROGRAM TERMINATION**********'); end
               break;
            end

            %update Lagragian multiplier
            [input]=update_lagran(itr,ip_cnt,input);

            %???
            [mesh]=kim_inversion2(input,mesh,fem);

            %plot updated model
            auto_contour_4d(itr,input,mesh,fem);
        end   %end for itr

        %FL debug
        if(input.debug==1); fprintf(input.debug_id,'\n===============\n'); end
        if(input.debug==1); fprintf(input.debug_id,'Out of iteration loop: objf=%f \n',fem.objf); end
        if(input.debug==1); fprintf(input.debug_id,'                       input.lagrn=%f \n',input.lagrn); end

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


