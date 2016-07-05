%
% ROUTINE DATA_2D: CREATE 2D ACQUISITION ARRAY
%
% Derived from function pushbutton1_Callback in original IP4DI data_2d.m file.
%
% Warning: does not enable borehole electrodes yet.
%          To use boreholes, a data file must be created manually and read in read_data.m
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: v1, October 15, 2015.

function [input]=data_2d(input)

disp(' ')
disp('=======================')
disp('=    ENTER DATA 2D    =')
disp('= create elect. array =')
disp('=======================')
disp(' ')


%% USER PARAMETERS
name_file=input.file_acqui_out
array_type=input.array_type;
if(array_type==1) disp('WENNER ARRAY'); end
if(array_type==2) disp('SCHLUMBERGER ARRAY'); end
if(array_type==3) disp('DIPOLE-DIPOLE ARRAY'); end
if(array_type==4) disp('POLE-DIPOLE ARRAY'); end
if(array_type==5) disp('POLE-POLE ARRAY'); end

number_electrodes=input.nb_electrodes
electrode_spacing=input.electrode_spacing
x_min=input.xelec_min
max_n=input.xN_max;
max_a=input.xA_max;
%%END USER PARAMETERS


% first create electrode coordinates
%FL: assuming a SURFACE configuration
%    (what about boreholes??)
xel=zeros(number_electrodes,1);
yel=zeros(number_electrodes,1);
for i=1:number_electrodes
    xel(i)=x_min+(i-1)*electrode_spacing;
    zel(i)=0;
end

pos_elec=[1:1:number_electrodes];

% determine all possible pairs of electrodes
max_pos_elec=nchoosek(pos_elec,2);
tmp_max_mes=size(max_pos_elec);


dd=[];
if array_type==1 %wenner

    tmp11=1;

    for i=1:tmp_max_mes(1)
        a=max_pos_elec(i,1);
        b=max_pos_elec(i,2);

        rest_electrodes=setdiff(pos_elec,a);
        rest_electrodes=setdiff(rest_electrodes,b);

        rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
        s_size=size(rest_electrodes);

        for j=1:s_size(1)
            m=rest_electrodes(j,1);
            n=rest_electrodes(j,2);

            if b-a<=max_a && m>a && n<b && m<n && abs(m-a)==abs(b-a) && abs(n-b)==abs(b-a)  
                dd(tmp11,1)=a;
                dd(tmp11,2)=b;
                dd(tmp11,3)=rest_electrodes(j,1);
                dd(tmp11,4)=rest_electrodes(j,2);

                tmp11=tmp11+1;
            end
        end
    end   


elseif array_type==2 %sch
% Schlumberger array
 
    tmp11=1;
    for i=1:tmp_max_mes(1)
        a=max_pos_elec(i,1);
        b=max_pos_elec(i,2);

        rest_electrodes=setdiff(pos_elec,a);
        rest_electrodes=setdiff(rest_electrodes,b);

        rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
        s_size=size(rest_electrodes);

        for j=1:s_size(1)
            m=rest_electrodes(j,1);
            n=rest_electrodes(j,2); 

            if b-a<=max_a && m>a && n<b && m<n && abs(m-a)<=max_n && abs(n-b)<=max_n  && abs(m-a)==abs(n-b)
                dd(tmp11,1)=a;
                dd(tmp11,2)=b;
                dd(tmp11,3)=rest_electrodes(j,1);
                dd(tmp11,4)=rest_electrodes(j,2);

                tmp11=tmp11+1;
            end
        end
    end       


elseif array_type==3 %dd
% Dipole-Dipole array

    tmp11=1;
    for i=1:tmp_max_mes(1)

        a=max_pos_elec(i,1);
        b=max_pos_elec(i,2);

        rest_electrodes=setdiff(pos_elec,a);
        rest_electrodes=setdiff(rest_electrodes,b);

        rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
        s_size=size(rest_electrodes);

        for j=1:s_size(1)
            m=rest_electrodes(j,1);
            n=rest_electrodes(j,2); 

            if b-a<=max_a && abs(min(m,n)-max(a,b))<=max_n && abs(m-n)==abs(b-a) && (min(m,n)>max(a,b))  

                dd(tmp11,1)=a;
                dd(tmp11,2)=b;
                dd(tmp11,3)=rest_electrodes(j,1);
                dd(tmp11,4)=rest_electrodes(j,2);

                tmp11=tmp11+1;

            end

        end
    end


elseif array_type==4 %pd
% Pole-Dipole array

    tmp11=1;
    for i=1:number_electrodes

        a=i;
        %b=max_pos_elec(i,2);

        rest_electrodes=setdiff(pos_elec,a);
        %rest_electrodes=setdiff(rest_electrodes,b);

        rest_electrodes=nchoosek(rest_electrodes,2); % Now this is M and N
        s_size=size(rest_electrodes);

        for j=1:s_size(1)
            m=rest_electrodes(j,1);
            n=rest_electrodes(j,2); 

            if abs(min(m,n)-a)<=max_n && abs(m-n)<=max_a && (max(m,n)<a || min(m,n)>a)

               dd(tmp11,1)=a;
               dd(tmp11,2)=0;
               dd(tmp11,3)=rest_electrodes(j,1);
               dd(tmp11,4)=rest_electrodes(j,2);

               tmp11=tmp11+1;

            end

        end

    end   %end for nb_electrodes


elseif array_type==5 %pp
% Pole-Pole array

    tmp11=1;
    for i=1:number_electrodes

        a=i;

        for j=1:number_electrodes

            if j~=i && abs(j-i)<=max_n
               dd(tmp11,1)=a;
               dd(tmp11,2)=0;
               dd(tmp11,3)=j;
               dd(tmp11,4)=0;

               tmp11=tmp11+1;
           end

        end

    end        

end   %end if array type


% Now add coordinates
clear data        
final=dd;

%Now assign coordinates in each electrode
max_num_mes=size(final)
max_pos_elec=final;   
 for i=1:max_num_mes(1)

     for j=1:number_electrodes

         if max_pos_elec(i,1)==pos_elec(j)
             data(i,1)=xel(j);
             data(i,2)=zel(j);
         end

         if max_pos_elec(i,2)==pos_elec(j)
             data(i,3)=xel(j);
             data(i,4)=zel(j);
         end       

         if max_pos_elec(i,3)==pos_elec(j)
             data(i,5)=xel(j);
             data(i,6)=zel(j);
         end 

         if max_pos_elec(i,4)==pos_elec(j)
             data(i,7)=xel(j);
             data(i,8)=zel(j);
         end        

     end

 end


% add virtual data
data(:,9)=10;
if (input.sip_flag==1 || input.ip_flag==1) data(:,10)=10; end

% write acqui and virtual data
%save(name_file,'data','-ascii');
fid=fopen(name_file,'wt');
  for i=1:max_num_mes(1)
      fprintf(fid,'%6.3f \t%6.3f \t%6.3f \t%6.3f \t%6.3f \t%6.3f \t%6.3f \t%6.3f \t%2.0f \t%2.0f \n',data(i,:));
  end
fclose(fid);

end   %end function
