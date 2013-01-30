    function final=info_save(stf,itr,ip_cnt,input,mesh,fem,final)
       
if stf==0
    final.num_param=mesh.num_param;
%     final.num_of_bor=mesh.num_of_bor;
    final.param_x=mesh.param_x;
    final.param_y=mesh.param_y;
    final.param_z=mesh.param_z;
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
    
    
    if input.bgr_res_flag~=0 &&input.time_lapse_flag==0 &itr>1  ;final.bgr_res_param=mesh.bgr_param;end

    if input.time_lapse_flag==1 ;final.num_files=input.num_files; final.d4_array_model_data=fem.d4_array_model_data; final.A=fem.A; end
    
    if itr>2 && input.acb_flag==1 
        final.ACB(:,itr-1)=fem.L1;
    end

    if itr>2  ;final.resolution_1=fem.resolution_matrix; end
    
    % Keep models
    final.res_param1(:,itr)=mesh.res_param2;
    final.acb(:,itr)=1;
    final.array_model_data=fem.array_model_data2;
    final.RMS(itr)=fem.rms_sum2;


    if input.ip_flag==1&& ip_cnt==1
        final.res_model=mesh.res_param2;
        final.res_final_data=final.array_model_data;
        final.all_res_model=final.res_param1;
        final.res_rms=final.RMS;
    end


    if input.ip_flag==1&& ip_cnt==2    
        final.chargeb=1000*(mesh.res_param2-mesh.res_final_param)./(mesh.res_param2);
        final.ip_model=final.res_param1;
        final.ip_rms(itr)=fem.rms_sum2;
    end

    if input.time_lapse_flag==1
        final.d4_res_param1(:,:,itr)=mesh.d4_res_param2;
    end

     




end

if stf==3
    
 if input.time_lapse_flag==0 %  I need on file for difference or backgroud inversion. All other files will be in mat form. 
    final.data(:,1)=input.ax;
    final.data(:,2)=input.ay;
    final.data(:,3)=input.az;
    final.data(:,4)=input.bx;
    final.data(:,5)=input.by;
    final.data(:,6)=input.bz;
    final.data(:,7)=input.mx;
    final.data(:,8)=input.my;
    final.data(:,9)=input.mz;
    final.data(:,10)=input.nx;
    final.data(:,11)=input.ny;
    final.data(:,12)=input.nz;
    final.data(:,13)=input.real_data;
    final.data(:,14)=fem.array_model_data;
     

    
    
%     write param and measurement	file
     
            datout=fopen('datout.inv','wt');
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if input.bgr_res_flag==0
                fprintf(datout,'INV_TYPE= %d DATAFILE= %s RMS= %f ITR= %d LGRN= %f NUM_PARAM= %d NUM_MEAS= %d\n',input.inv_flag,input.mes_in,fem.rms_sum2,itr-1,input.lagrn,mesh.num_param,input.num_mes);
                if input.dc_flag==1
                    for i=1:mesh.num_param
                        fprintf(datout,'%f %f %f %f\n',mesh.param_x(i),mesh.param_y(i),-mesh.param_z(i),mesh.res_param1(i));
                    end 
                    fprintf(datout,'---------------------------------\n');

                        for i=1:input.num_mes
                        fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' ...
                        ,input.m_array_type(i), final.data(i,1), final.data(i,2), final.data(i,3),...
                         final.data(i,4), final.data(i,5), final.data(i,6),...
                          final.data(i,7), final.data(i,8), final.data(i,9),...
                           final.data(i,10), final.data(i,11), final.data(i,12),...
                            final.data(i,13),final.data(i,14));
                        end
                elseif input.ip_flag==1
                    % Not yet supported
                elseif input.sip_flag==1
                    for i=1:mesh.num_param
                        fprintf(datout,'%f %f %f %f %f\n',...
                            mesh.param_x(i),mesh.param_y(i),-mesh.param_z(i),...
                            abs(mesh.res_param1(i)),1000*(  atan2(imag(mesh.res_param1(i)),real(mesh.res_param1(i))   ) ) );
                    end 
                     fprintf(datout,'---------------------------------\n'); 

                    for i=1:input.num_mes
                        fprintf(datout,'%d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n' ...
                        ,input.m_array_type(i), final.data(i,1), final.data(i,2), final.data(i,3),...
                         final.data(i,4), final.data(i,5), final.data(i,6),...
                          final.data(i,7), final.data(i,8), final.data(i,9),...
                           final.data(i,10), final.data(i,11), final.data(i,12),...
                            real(final.data(i,13)),imag(final.data(i,13)),...
                            real(final.data(i,14)),imag(final.data(i,14)));
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fclose(datout);
        end

    
end

    end