% 	  /* ******************************************************************
% 	   *                             RMS( )                               *
% 	   *------------------------------------------------------------------*
% 	   * This function finds the RMS error and ends the procedure         *
% 	   ******************************************************************** */


function [ex,fem,mesh]=rms_4d(itr,ip_cnt,input,mesh,fem)

global textt
fem.rms_sum1=0;
ex=0;
tsum=0;

% Create real_data in one 
tmp_d4=[];
for i=1:input.num_files;
    tmp_d4=[tmp_d4;input.d4_real_data(:,i)];
end
        

% Find RMS error
for i=1:input.num_mes*input.num_files
    tsum=tsum+(fem.e(i)^2)./(tmp_d4(i)^2);
    
end
fem.rms_sum1=sqrt(tsum/input.num_mes*input.num_mes)*100;
disp(sprintf('RMS=>%f',fem.rms_sum1)); 


% Update for first iteration
if (itr==1)
    fem.rms_sum2=fem.rms_sum1;
end



if (fem.rms_sum1>fem.rms_sum2 && itr>1)
    ka=sprintf('Divergence occured. Old RMS=>%f New Rms=>%f',fem.rms_sum2,fem.rms_sum1);
    set(textt,'String',ka);drawnow;
    % If we have IP data update everything to continue with IP inversion 
    if input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==2)
        ex=1;
    elseif input.ip_flag==1 && ip_cnt==1
        ex=2;        
    end
        
end


tmp11=abs((fem.rms_sum2-fem.rms_sum1)/fem.rms_sum2);

if ( tmp11<(input.conv_rate/100)  && itr>1)
    set(textt,'String','Very slow convergence');drawnow
    % If we have IP data update everything to continue with IP inversion
    if input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==2)
        ex=1;
    elseif input.ip_flag==1 && ip_cnt==1
        ex=2;     
    end
end

if ex==0
    mesh.d4_res_param2=mesh.d4_res_param1;
    fem.rms_sum2=fem.rms_sum1;
    fem.d4_array_model_data2=fem.d4_array_model_data;
end


if itr==input.itn+1 && ex==0
            set(textt,'String','********PROGRAM TERMINATION**********');drawnow
            
            
            % If we have IP data update everything to continue with IP inversion
            if input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==2)
             ex=1;
            elseif input.ip_flag==1 && ip_cnt==1
             ex=2;            
            end
end





end
