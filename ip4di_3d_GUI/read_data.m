function [input]=read_data(input)




tmp_d=importdata(input.mes_in);


% This is the common part
input.ax=tmp_d(:,1);
input.ay=tmp_d(:,2);
input.az=tmp_d(:,3);
input.bx=tmp_d(:,4);
input.by=tmp_d(:,5);
input.bz=tmp_d(:,6);
input.mx=tmp_d(:,7);
input.my=tmp_d(:,8);
input.mz=tmp_d(:,9);
input.nx=tmp_d(:,10);
input.ny=tmp_d(:,11);
input.nz=tmp_d(:,12);
    



input.real_data=tmp_d(:,13);


input.num_mes=length(input.ax);



%%%%%%%%%%%%%%%%%%%%%%%%%
%%%IP DATA or SIP%%%%%%
try
    input.ip_data=tmp_d(:,14); 
    
    
    if input.sip_flag==1
        disp('SIP data found');
        input.ip_num=1;
        input.real_data=complex(input.real_data,input.ip_data);
        % if are in mag and phase
%         real_part=input.real_data.*(cos(input.ip_data./1000));
%         imag_part=input.real_data.*(sin(input.ip_data./1000));
%         input.real_data=complex(real_part,imag_part);
         
    elseif input.ip_flag==1
        disp('IP data found');
        input.ip_num=2;
        input.ip_data=input.ip_data./1000;% Loke has it in msec
    else
        input.stdev_error=input.ip_data; % Currently on fow DC
	input.ip_num=1;
    end
catch
    disp('No IP or SIP data found');
    input.ip_num=1;
    input.stdev_error=ones(input.num_mes,1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here check for IP or SIP and STDEV
if input.sip_flag==1 || input.ip_flag==1
    try
        input.stdev_error=tmp_d(:,15); 
    catch
        input.stdev_error=ones(input.num_mes,1);
    end
end


clear tmp_d;

input.Wd=diag(1./input.stdev_error).^2;




        

% Try to find number of electrodes used




%--------------------------------------------------------------------------
    % search for how many electrodes
    error=0;
    b_el=0;
    n_el=0;
    input.array_type=1;
    ax=length(unique(input.ax));
    ay=length(unique(input.az));
    az=length(unique(input.az));
    
    bx=length(unique(input.bx));
    by=length(unique(input.bz));
    bz=length(unique(input.bz));

    
    mx=length(unique(input.mx));
    my=length(unique(input.mz));
    mz=length(unique(input.mz));

    
    nx=length(unique(input.nx));
    ny=length(unique(input.nz));
    nz=length(unique(input.nz));
    
    if ax==1 && ay==1 && az==1
        disp('Error. A electrode must be there ALWAYS');
        error=1;
    end
    
    if mx==1 && my==1 && mz==1
        disp('Error. M electrode must be there ALWAYS');
        error=1;
    end
    
    if bx==1 && by==1 && bz==1
        b_el=1;
        error=0;
    end
    
    if nx==1 && ny==1 && nz==1
        n_el=1;        
        error=0;
    end
    
    
   if b_el==1 && n_el==0
       disp('Pole Dipole array?');
       input.elec_array_type=2;
   elseif b_el==1 && n_el==1
       disp('Pole Pole array?');
       input.elec_array_type=3;
   elseif b_el==0 && n_el==1
       disp('Not tested array');
       pause
   else
       disp('4 electrode array?');
       input.elec_array_type=1;
   end



if error==1
    pause
end

      
input.m_array_type(1:input.num_mes,1)=input.elec_array_type;



end