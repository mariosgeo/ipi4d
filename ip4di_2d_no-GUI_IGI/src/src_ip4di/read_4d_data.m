function [input]=read_4d_data(input)

 % read data file
 tmp=importdata(input.mes_in);
 input.num_files=length(tmp);
 disp(sprintf('NB OF TIME-LAPSE EXPERIMENTS = %i',input.num_files))

 try
    % recast data in a 3rd-order tensor [i_row=i_data,i_col=elec/data,i_file]
    for i=1:input.num_files
        %tmp2(:,:,i)=importdata(tmp{i});
        tmp2(:,:,i)=tmp{i};
    end
 catch
    disp('this version does not support files with different number of measurements on each file. Future may...');
    return
 end   %end try

 input.ax=tmp2(:,1,1);
 input.az=abs(tmp2(:,2,1));
 input.bx=tmp2(:,3,1);
 input.bz=abs(tmp2(:,4,1));
 input.mx=tmp2(:,5,1);
 input.mz=abs(tmp2(:,6,1));
 input.nx=tmp2(:,7,1);
 input.nz=abs(tmp2(:,8,1));

 % nb of data for one time-lapse experiment
 input.num_mes=length(input.ax);
 disp(sprintf('NB OF MEASUREMENTS = %i',input.num_mes))
 disp(' ')

 % first dataset (real part)
 real_data=tmp2(:,9,1);


 %%%%%%%%%%%%%%%%%%%%%%%
 %   IP DATA or SIP    %
 try %if in first file, there is a 10th column, then check for IP or SIP
    ip_data=tmp2(:,10,1); 

    if input.sip_flag==1
        disp('SIP data found');
        input.ip_num=1;
            for i=1:input.num_files
                ip_data=tmp2(:,10,i); 
                real_data=tmp2(:,9,i);

                % time-lapse data
                input.d4_real_data(:,i)=complex(real_data,ip_data);     

                % if are in mag and phase
                % real_part=real_data(:,i).*(cos(ip_data(:,i)./1000));
                % imag_part=real_data(:,i).*(sin(ip_data(:,i)./1000));
                % input.real_data(:,i)=complex(real_part,imag_part);
            end
         
    elseif input.ip_flag==1
        disp('IP data found');
        input.ip_num=2;
            for i=1:input.num_files        
                input.d4_ip_data(:,i)=tmp2(:,10,i)./1000; 
                % time-lapse data
                input.d4_real_data(:,i)=complex(real_data,ip_data);      
            end
    else
    % if not SIP or IP, then the 10th col. corresponds to uncertainties
        
        for i=1:input.num_files
            % time-lapse data
            input.d4_real_data(:,i)=tmp2(:,9,i);     
        end    
        
        for i=1:input.num_files
            % uncertainties
            input.d4_stdev_error(:,i)=tmp2(:,10,i);   %MK: Currently on fow DC
        end
        input.ip_num=1;
    end

catch
% no 10th col. found
    disp('No IP or SIP data found');
    input.ip_num=1;
    for i=1:input.num_files
        % time-lapse data
        input.d4_real_data(:,i)=tmp2(:,9,i);     
    end  
    
    input.d4_stdev_error=ones(input.num_mes,input.num_files);   %MK: I need to change it to take accound the 4D
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Here check for IP or SIP uncertainties (11th col.)
if input.sip_flag==1 || input.ip_flag==1
    try
        for i=1:input.num_files
            input.d4_stdev_error(:,i)=tmp2(:,11,i);
        end        
    catch
        %FL: changed stdev_error to d4_stdev_error (correct, I think...)
        input.d4_stdev_error=ones(input.num_mes,input.num_files);
    end
end


% define single stdev_error (FL?)
input.stdev_error=ones(input.num_mes,1);

%input.Wd=diag(input.stden_error); 
input.Wd_d4=sparse(eye(input.num_files*input.num_mes)).^2;



%--------------------------------------------------------------------------
% Try to find number of electrodes used

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
       disp('Pole Dipole array');
       input.elec_array_type=2;
   elseif b_el==1 && n_el==1
       disp('Pole Pole array');
       input.elec_array_type=3;
   elseif b_el==0 && n_el==1
       disp('Not tested array');
       pause
   else
       disp('4 electrode array');
       input.elec_array_type=1;
   end

 if error==1
    pause
 end

 input.m_array_type(1:input.num_mes,1)=input.elec_array_type;

end

