% 	  /* *******************************************************************
% 	   *                      PROP_CHANGE()                                *
% 	   *-------------------------------------------------------------------*
% 	   * This function changes the property of the elements according to   *
% 	   * the outcome of the iteration in order to have updated property    *
% 	   * for the next iteration.                                           *
% 	   *                                                                   *
% 	   ********************************************************************* */



function [fem,mesh,input]=prop_change(itr,input,mesh,fem)


if input.bgr_res_flag~=0%&& itr>1
    
     disp ('              --------- READ BACKGROUND MODEL --------');
      
    % Now we have two options, DC or IP or SIP
    
    [misc,tin,misc,misc,misc,tmp15,misc,nitr,misc,tmp16_lag,misc,tpar,misc,tmes]=...
                                    textread(input.par_nm,'%s %d %s %s %s %f %s %d %s %f %s %d %s %d',1);
    
    if tpar~=mesh.num_param
        disp('Background Model not Compatible with these Parameters !!!') 
        disp('PROGRAM TERMINATION !!!!')
        pause
    end
    
    if tmes~=input.num_mes
        disp('Background measurments not Compatible with these measurments !!!') 
        disp('PROGRAM TERMINATION !!!!')
        pause
    end
    
       
    tmp_d=importdata(input.par_nm,' ' ,1);
%     tmp_d.data(tpar,:)=0;
    x_center=tmp_d.data(1:tpar,1);
    y_center=tmp_d.data(1:tpar,2);
    mesh.bgr_param=tmp_d.data(1:tpar,3);
   

    if input.ip_flag==1 
         tmp=tmp_d.data(1:tpar,4); % this is chargebility         
    end
    
    if input.sip_flag==1
        tmp=tmp_d.data(1:tpar,4); % this is in mrads
        real_part=mesh.bgr_param.*(cos(tmp./1000));
        imag_part=mesh.bgr_param.*(sin(tmp./1000));
        mesh.bgr_param=complex(real_part,imag_part);
    end
    


         %This is the background array data
         %[bgr_m_array_type,tax,taz,tbx,tbz,tmx,tmz,tnx,tnz,real_data0,array_model_data_bgr]=textread(par_nm,'%d %f %f %f %f %f %f %f %f %f %f','headerlines',tpar+1);
         tmp_d=importdata(input.par_nm,' ' ,tpar+1);
         bgr_m_array_type=tmp_d.data(:,1); 
         tax=tmp_d.data(:,2);  
         taz=tmp_d.data(:,3);  
         tbx=tmp_d.data(:,4);  
         tbz=tmp_d.data(:,5);  
         tmx=tmp_d.data(:,6);  
         tmz=tmp_d.data(:,7);  
         tnx=tmp_d.data(:,8);  
         tnz=tmp_d.data(:,9);  
         
         if input.dc_flag==1
             input.real_data0=tmp_d.data(:,10);  
             input.array_model_data_bgr=tmp_d.data(:,11); 
         elseif input.ip_flag==1
             input.real_data0=tmp_d.data(:,10);  
             input.ip_data0=tmp_d.data(:,11);
             input.array_model_data_bgr=tmp_d.data(:,12); 
             input.ip_array_model_data_bgr=tmp_d.data(:,13);             
         elseif input.sip_flag==1
             input.real_data0=complex(tmp_d.data(:,10),tmp_d.data(:,11));
             input.array_model_data_bgr=complex(tmp_d.data(:,12),tmp_d.data(:,13));
         end


    
   mesh.res_param1=mesh.bgr_param;
   mesh.res_param2=mesh.res_param1;
   
  
    
end


% 
% for i=1:mesh.num_param
% 
%     for j=1:mesh.num_elements
% 
%         if(mesh.icon(4,j)==abs(i))  mesh.prop(j)=mesh.res_param1(i); end
%     end
% end


    for i=1:mesh.num_param
        ind= mesh.icon(4,:)==i;
        mesh.prop(ind)=mesh.res_param1(i);

    end
end
