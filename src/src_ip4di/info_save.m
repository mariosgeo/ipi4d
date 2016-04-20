%
% Function INFO_SAVE: save variables of interest for inversion
%                     - at current iteration (update from one iteration to another)
%                     - along iterations (store variables vs. it.)
%
% Author: Marios Karaoulis, Colorado School of Mines

% Modified by Francois Lavoue', Colorado School of Mines, September 2015.
% Modifications include:
% 1) some comments for more readability,
% 2) saving intermediate data and models vs. it.
% 3) change reference to 'mesh.res_param2' (model at previous iteration) in 'mesh.res_param1' (current model)
%    (except for IP and time-lapse inversion)
%    and idem for data and RMS (given that, when info_save is called, RMS, data and model consistently refer
%    to current iteration, before model update).
% 4) some change of notations.
%    - 'final.res_param1' is now 'final.res_param1_vs_it'

function final=info_save(stf,itr,ip_cnt,input,mesh,fem,final)

if stf==0
% Intermediate iteration case:
% store data and parameters in variables

    %
    % SAVE CURRENT VARIABLES
    %
    final.num_param=mesh.num_param;
    final.num_of_bor=mesh.num_of_bor;
    final.param_x=mesh.param_x;
    final.param_y=mesh.param_y;
    % final.borehole_no_surface_flag=input.borehole_no_surface_flag;
    final.itn=itr;
    final.lagrn=input.lagrn;
    final.inv_flag=input.inv_flag;
    final.bgr_res_flag=input.bgr_res_flag;
    % final.par_nm=input.par_nm;
    final.acb_flag=input.acb_flag;
    final.mes_in=input.mes_in;
    final.num_mes=input.num_mes;
    final.time_lapse_flag=input.time_lapse_flag;
    final.ip_flag=input.ip_flag;
    final.sip_flag=input.sip_flag;
    final.dc_flag=input.dc_flag;

    % update current synthetic data
    final.array_model_data=fem.array_model_data;

    %in case of difference inversion, save background model
    if input.bgr_res_flag~=0 && input.time_lapse_flag==0  && itr>1; final.bgr_res_param=mesh.bgr_param; end


    %
    % SAVE VARIABLES VS ITERATIONS
    %
    final.RMS_vs_it(itr,:)=[itr fem.rms_crit ...
                                fem.objf fem.objf_data fem.objf_model input.lagrn*fem.objf_model ...
                                fem.rms fem.nrms fem.wrms fem.mrms];   %save misfits vs. it
    final.res_param1_vs_it(:,itr)=mesh.res_param1;        %save models vs. it.
    final.model_data_vs_it(:,itr)=fem.array_model_data;   %save synthetic data vs. it
    final.acb(:,itr)=1;   %FL: ??

    %in case of ACB, save Lagrangian distributions vs. it.
    if input.acb_flag==1
        final.ACB(:,itr)=fem.L1;
    end

    %in case of time-lapse inversion,
    if input.time_lapse_flag==1
       final.num_files=input.num_files;                  %save nb of data files
       final.d4_res_param1_vs_it(:,:,itr)=mesh.d4_res_param1;        %store all models vs. it.
       final.d4_model_data_vs_it(:,:,itr)=fem.d4_array_model_data;   %save synthetic data vs. it
    end


    % IP inversion: 1st component
    if input.ip_flag==1 && ip_cnt==1
        final.res_model=mesh.res_param;
        final.res_final_data=final.array_model_data;
        final.all_res_model=mesh.res_param;
        final.res_rms=final.RMS_vs_it(:,2);
    end

    % IP inversion: 2nd component
    if input.ip_flag==1 && ip_cnt==2
        final.chargeb=1000*(mesh.res_param2-mesh.res_final_param)./(mesh.res_param2);
        final.ip_model=mesh.res_param2;
        final.ip_rms(itr)=fem.rms_crit2;
    end


    % Print intermediate infos (misfit values vs. it.)
    if itr==1
    % open info file for first time
       dattxt=fopen(input.file_info,'wt');
       if input.time_lapse_flag==0
          fprintf(dattxt,'NB PARAM = %d, NB DATA = %d, INV TYPE = %d (%s), ACB = %d \n\n', ...
                          mesh.num_param,input.num_mes,input.inv_flag,input.inv_name,input.acb_flag);
       else   %time-lapse output
          fprintf(dattxt,'NB PARAM = %d, NB DATA = %d, INV TYPE = %d (%s), ACB = %d, GAMMA = %f \n\n', ...
                          mesh.num_param,input.num_mes,input.inv_flag,input.inv_name,input.acb_flag,input.gamma);
       end

       % print data file
       fprintf(dattxt,'DATA FILE = %s\n\n', input.mes_in);

       if input.image_guidance>0
          fprintf(dattxt,'TI FILE = %s\n\n', input.training_image)
       end

       if input.time_lapse_flag==0
          fprintf(dattxt,'ITR \t C(m) \t C_D(m) \t C_M(m) \t l*C_M(m) \t l*C_M/C_D \t nRMS (%%) \t LAMBDA \n');
       else   %time-lapse output (FL: syntax is not great, but output fits in a 13" screen...)
          %fprintf(dattxt,'ITR \t\b\b\b C(m) \t\b  C_D(m) \t\b\b\b\b\b C_M(m) \t\b l*C_M(m) \t\b\b\b\b\b l*C_M/C_D \t\b\b\b C_T(m) \t\b\b\b g*C_T(m) \t\b\b\b g*C_T/C_D \t\b\b\b nRMS (%%) \t\b\b\b LAMBDA \n');
          %fprintf(dattxt,'ITR    C(m)        C_D(m)      C_M(m)       l*C_M(m)     l*C_M/C_D    C_T(m)     g*C_T(m)   g*C_T/C_D  nRMS (%)   LAMBDA\n');
          fprintf(dattxt,'ITR\t C(m)\t\t C_D(m)\t\t C_M(m)\t\t l*C_M(m)\t l*C_M/C_D\t C_T(m)\t\t g*C_T(m)\t g*C_T/C_D\t nRMS (%%)\t LAMBDA \n');

       end

    else
    % add new info line to file
       dattxt=fopen(input.file_info,'a');
    end

    % print info for current iteration
    %fprintf(dattxt,'INV_TYPE= %d DATAFILE= %s RMS= %f ITR= %d LGRN= %f NUM_PARAM= %d NUM_MEAS= %d\n',...
    %               input.inv_flag,input.mes_in,fem.nrms2,itr-1,input.lagrn,mesh.num_param,input.num_mes);
    if input.time_lapse_flag==0
       %               ITR  C(m)  C_D   C_M l*C_M l*C_M/C_D nRMS LGRN
       fprintf(dattxt,'%i \t %f \t %f \t %f \t %f \t %f \t %f \t %f \n', ...
                       itr-1,fem.objf,fem.objf_data,...
                       fem.objf_model,input.lagrn*fem.objf_model,input.lagrn*fem.objf_model/fem.objf_data,...
                       fem.nrms,input.lagrn);

    else   %time-lapse output
       %               ITR   C(m) C_D  C_M  l*C_M l*C_M/C_D C_T g*C_T g*C_T/C_D  nRMS  LAMBDA
       %fprintf(dattxt,'%i \t\b\b\b %f \t\b  %f \t\b\b\b\b\b %f \t\b %f \t\b\b\b\b\b %f \t\b\b\b %f \t\b\b\b %f \t\b\b\b %f \t\b\b\b %f \t\b\b\b %f \n', ...
       %fprintf(dattxt,'%i3    %f    %f    %f    %f    %f    %f    %f    %f    %f    %f\n', ...
       fprintf(dattxt,'%i\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f\t %f \n', ...
                       itr-1, fem.objf, fem.objf_data,...
                       fem.objf_model, input.lagrn*fem.objf_model, input.lagrn*fem.objf_model/fem.objf_data,...
                       fem.objf_tl,    input.gamma*fem.objf_tl,    input.gamma*fem.objf_tl/fem.objf_data,...
                       fem.nrms, input.lagrn);
    end

    %close info file
    fclose(dattxt);

end   %end if stf==0


if stf==1
% Intermediate iteration case:
% store resolution and Jacobian matrices
% (after 'invert_cntr.m where these amtrices are computed)

    %in case of non-time-lapse inversion,
    if input.time_lapse_flag==0
       final.resolution_vs_it(:,itr)=diag(fem.resolution_matrix);             %save resolution matrix vs. it. (diagonal only)
       final.resolution_data_vs_it(:,itr)=diag(fem.resolution_matrix_data);   %save data resolution matrix vs. it. (diagonal only)
       final.jacobian_vs_it(:,:,itr)=fem.array_jacobian;                      %save Jacobian matrix vs. it.
    %FL: and what about time-lapse inversion?
    end

end


if stf==0 && input.print_intermediate_results==1
% Intermediate iteration case:
% print data and parameters in files

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   if input.time_lapse_flag==0 
   % NON-TIME-LAPSE INVERSION

      %FL: Data and models are stored in column format,
      %    so we have to rewrite each time data and models for all iterations
      %    This is done with the 'save' command instead of fopen/fprintf/fclose.
      %datout=fopen(input.file_data_inv_inter,'wt');
      %modout=fopen(input.file_model_inter,'wt');

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if input.bgr_res_flag==0
      %case where we don't use a background model (FL: and what ifwe do?...)

         if input.dc_flag==1
         % DC data

            % output all intermediate models until current iteration
            Mts=[mesh.param_x' mesh.param_y' final.res_param1_vs_it];
            save(input.file_model_inter,'Mts','-ascii');

            % output obs vs synthetic data for all iterations
            Mts=[ones(input.num_mes,1),input.ax,input.az,input.bx,input.bz,input.mx,input.mz,input.nx,input.nz,...
                 input.real_data,final.model_data_vs_it];
            save(input.file_data_inv_inter,'Mts','-ascii')

         elseif input.ip_flag==1
         % IP data (not yet supported)

         elseif input.sip_flag==1
         % SIP data and models

            if input.cmplx_format==1
            %save real and imaginary part of reconstructed resistivity at all iterations
               Mts=[mesh.param_x,mesh.param_y];
               for iitr=1:itr
                   Mts=[Mts , real(final.res_param1_vs_it(:,itr)) , imag(final.res_param1_vs_it(:,itr))];
               end
            elseif input.cmplx_format==2
            %save amplitude and phase of reconstructed resistivity at all iterations
               Mts=[mesh.param_x,mesh.param_y];
               for iitr=1:itr
                   Mts=[Mts , abs(mesh.res_param1_vs_it(:,itr)) , 1000*atan2(imag(mesh.res_param1_vs_it(:,itr)),real(mesh.res_param1_vs_it(:,itr))) ];
               end
            end
            save(input.file_model_inter,'Mts','-ascii');

            % output obs vs synthetic data for all iterations
            if input.cmplx_format==1
            %save real and imaginary part of synthetic apparent resistivity
               Mts=[ ones(input.num_mes,1),input.ax,input.az,input.bx,input.bz,input.mx,input.mz,input.nx,input.nz,...
                     real(input.real_data),imag(input.real_data),real(final.model_data_vs_it),imag(final.model_data_vs_it) ];
            elseif input.cmplx_format==2
            %save amplitude and phase of synthetic apparent resistivity
               Mts=[ ones(input.num_mes,1),input.ax,input.az,input.bx,input.bz,input.mx,input.mz,input.nx,input.nz,...
                     abs(       input.real_data),1000*atan2(imag(       input.real_data),real(       input.real_data)),...
                     abs(final.model_data_vs_it),1000*atan2(imag(final.model_data_vs_it),real(final.model_data_vs_it)) ];
            end
            save(input.file_data_inv_inter,'Mts','-ascii')
         end   %end if data_flag
      end   %end if bgr_res_flag==0


   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   elseif input.time_lapse_flag==1
   % TIME-LAPSE INVERSION

   %FL: TO DO

   end   %end non-time-lapse vs time-lapse inversion

end



if stf==3   %Final iteration case:
% write parameter and measurement files

   if input.time_lapse_flag==0
   % non-time-lapse inversion 
   %MK: I need on file for difference or backgroud inversion. All other files will be in mat form (FL: ??)

      datout=fopen(input.file_data_inv_out,'wt');
      modout=fopen(input.file_model_out,'wt');
      modint=fopen(input.file_model_interp,'wt');

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      if input.bgr_res_flag==0

         if input.dc_flag==1
         % DC data

            % output final model
            for i=1:mesh.num_param
                fprintf(modout,'%f %f %f\n',mesh.param_x(i),mesh.param_y(i),mesh.res_param1(i));
            end
            % output obs vs synthetic data
            for i=1:input.num_mes
                fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f\n',...
                               1,input.ax(i),input.az(i),input.bx(i),input.bz(i),input.mx(i),input.mz(i),input.nx(i),input.nz(i),...
                               input.real_data(i),fem.array_model_data(i));
            end

         elseif input.ip_flag==1
         % IP data (not yet supported)

         elseif input.sip_flag==1
         % SIP data and models

            % output final model
            for i=1:mesh.num_param
                if input.cmplx_format==1
                %save real and imaginary part of reconstructed resistivity
                   fprintf(modout,'%f %f %f %f\n',...
                           mesh.param_x(i),mesh.param_y(i),real(mesh.res_param1(i)),imag(mesh.res_param1(i)) );
                elseif input.cmplx_format==2
                %save amplitude and phase of reconstructed resistivity
                   fprintf(modout,'%f %f %f %f\n',...
                           mesh.param_x(i),mesh.param_y(i),abs(mesh.res_param1(i)),1000*atan2(imag(mesh.res_param1(i)),real(mesh.res_param1(i))) );
                end
            end

            % eventually, output image-guided interpolation of final model
            if input.image_guided_interpolation==1
               % perform image-guided interpolation of final model
               mesh.res_param_interp=image_guided_interpolation(input,mesh,mesh.res_param1);

               % output
               for i=1:mesh.num_param
                   if input.cmplx_format==1
                   %save real and imaginary part of reconstructed resistivity
                      fprintf(modint,'%f %f %f %f\n',...
                              mesh.param_x(i),mesh.param_y(i),real(mesh.res_param_interp(i)),imag(mesh.res_param_interp(i)) );
                   elseif input.cmplx_format==2
                   %save amplitude and phase of reconstructed resistivity
                      fprintf(modint,'%f %f %f %f\n',...
                              mesh.param_x(i),mesh.param_y(i),abs(mesh.res_param_interp(i)),1000*atan2(imag(mesh.res_param_interp(i)),real(mesh.res_param_interp(i))) );
                   end
               end
            end   %end if image-guided interpolation

            % output obs vs synthetic data
            for i=1:input.num_mes
                if input.cmplx_format==1
                %save real and imaginary part of synthetic apparent resistivity
                   fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f %f %f\n',...
                                  1,input.ax(i),input.az(i),input.bx(i),input.bz(i),input.mx(i),input.mz(i),input.nx(i),input.nz(i),...
                                  real(input.real_data(i)),imag(input.real_data(i)),real(fem.array_model_data(i)),imag(fem.array_model_data(i)) );
                elseif input.cmplx_format==2
                %save amplitude and phase of synthetic apparent resistivity
                   fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f %f %f\n',...
                                  1,input.ax(i),input.az(i),input.bx(i),input.bz(i),input.mx(i),input.mz(i),input.nx(i),input.nz(i),...
                                  abs(     input.real_data(i)),1000*atan2(imag(     input.real_data(i)),real(     input.real_data(i))),...
                                  abs(fem.array_model_data(i)),1000*atan2(imag(fem.array_model_data(i)),real(fem.array_model_data(i))) );
                end
            end

         end   %end if data_flag
      end   %end if bgr_res_flag==0
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

      fclose(datout);
      fclose(modout);
      fclose(modint);
   end   %end if time_lapse_flag==0 (FL: and what about time-lapse results??)
          
end   %end if stf=3 (=end of inversion)


end   %end function info_save


