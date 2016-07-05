function [mesh]=create_mesh3(input,manual_flag)




mesh=find_probe_spacing(input);
if manual_flag==0;
    mesh=find_max_n(input,mesh);
else
   mesh.max_n=length(manual_flag);
   mesh.depth_n=manual_flag;
end
mesh=calc_param(input,mesh);
mesh=create_mesh_gen(input,mesh);
% calculate_half_band
mesh=create_mesh_grid(input,mesh);
mesh=geometrical_factor(input,mesh);


mesh.init_res_param(1:mesh.num_param)=mesh.mean_res;
mesh.res_param1(1:mesh.num_param,1)=mesh.mean_res;
mesh.res_param2(1:mesh.num_param,1)=mesh.mean_res;
mesh.prop(1:mesh.num_elements)=mesh.mean_res;

if input.time_lapse_flag==1
    mesh.d4_res_param1=mesh.mean_res*ones(mesh.num_param,input.num_files);
    mesh.d4_res_param2=mesh.mean_res*ones(mesh.num_param,input.num_files);   
end

end




function mesh=geometrical_factor(input,mesh)
% in

% include first topography in ax and az coordinates
% for i=1:input.num_mes
%     for j=1:length(mesh.anaglyfo_data(:,1))
%         if input.ax(i)==mesh.anaglyfo_data(j,1)
%             input.az(i)=input.az(i)-anaglyfo_data(j,2);
%         end
%         if input.bx(i)==mesh.anaglyfo_data(j,1)
%             input.bz(i)=input.bz(i)-anaglyfo_data(j,2);
%         end
%         if input.mx(i)==mesh.anaglyfo_data(j,1)
%             input.mz(i)=input.mz(i)-anaglyfo_data(j,2);
%         end
%         if input.nx(i)==mesh.anaglyfo_data(j,1)
%             input.nz(i)=input.nz(i)-anaglyfo_data(j,2);
%         end
%     end
% end


%Find the geometrical factor
VAM=zeros(input.num_mes,1);
VAN=zeros(input.num_mes,1);
VBM=zeros(input.num_mes,1);
VBN=zeros(input.num_mes,1);




for i=1:input.num_mes
    if(input.m_array_type(i)==1 || input.m_array_type(i)==2)
         aflag=1; bflag=1; mflag=1; nflag=1;
    end
             %/* flag for Pole-dipole array */
    if(input.m_array_type(i)==2)
        aflag=1; bflag=0; mflag=1; nflag=1;
    end
		 %/* flag for Pole-pole */
    if(input.m_array_type(i)==3)
          aflag=1; mflag=1; bflag=0; nflag=0;
    end
    
    
if(aflag*mflag~=0)
    r1=sqrt( (input.ax(i)-input.mx(i))*(input.ax(i)-input.mx(i)) + (input.az(i)-input.mz(i))*(input.az(i)-input.mz(i)) );
    r2=sqrt( (input.ax(i)-input.mx(i))*(input.ax(i)-input.mx(i)) + (-input.az(i)-input.mz(i))*(-input.az(i)-input.mz(i)) );
    VAM(i)=(r1+r2)./(r1.*r2.*4*pi);
end

if(aflag*nflag~=0)
    r1=sqrt( (input.ax(i)-input.nx(i))*(input.ax(i)-input.nx(i)) + (input.az(i)-input.nz(i))*(input.az(i)-input.nz(i)) );
    r2=sqrt( (input.ax(i)-input.nx(i))*(input.ax(i)-input.nx(i)) + (-input.az(i)-input.nz(i))*(-input.az(i)-input.nz(i)) );
    VAN(i)=(r1+r2)./(r1.*r2.*4*pi); 
end

if(bflag*mflag~=0)
    r1=sqrt( (input.bx(i)-input.mx(i))*(input.bx(i)-input.mx(i)) + (input.bz(i)-input.mz(i))*(input.bz(i)-input.mz(i)) );
    r2=sqrt( (input.bx(i)-input.mx(i))*(input.bx(i)-input.mx(i)) + (-input.bz(i)-input.mz(i))*(-input.bz(i)-input.mz(i)) );
    VBM(i)= (r1+r2)./(r1.*r2.*4*pi); 
end

if(bflag*nflag~=0)
    r1=sqrt( (input.bx(i)-input.nx(i))*(input.bx(i)-input.nx(i)) + (input.bz(i)-input.nz(i))*(input.bz(i)-input.nz(i)) );
    r2=sqrt( (input.bx(i)-input.nx(i))*(input.bx(i)-input.nx(i)) + (-input.bz(i)-input.nz(i))*(-input.bz(i)-input.nz(i)) );
    VBN(i)= (r1+r2)./(r1.*r2.*4*pi); 
end
end

mesh.geofac=1./(VAM-VBM-VAN+VBN);
clear VAM;clear VBM;clear VAN;clear VBN;clear r1;clear r2;



if input.data_type==2 
input.real_data=input.real_data.*mesh.geofac;
%pause
end






if input.time_lapse_flag==1
    ind=find (input.d4_real_data<0);
    if size(ind,1)~=0 && size(ind,2)~=-0
        display('Negative resistivitys in measurments');
        ind
        pause
    end   
    mesh.mean_res=mean(mean(10.^(log10(input.d4_real_data))));  
    
elseif input.time_lapse_flag==0
    ind=find (input.real_data<0);
    if size(ind,1)~=0 && size(ind,2)~=-0
        display('Negative resistivitys in measurments');
        ind
        pause
    end    
    mesh.mean_res=(mean(10.^(log10(input.real_data))));
end

end



% The main idea for this function is to calculate the basic dimensions of
% the mesh. I need to find the surface electrodes, the borehole electrodes,
% the positions of the electrodes etc. If for some reason a borehole is on
% edge or i hane no surface electrodes i must ensure to add the neccessaty
% points.

function mesh=find_probe_spacing(input)





%%%%%%%%%%%%%%%%%%%%%%PRE-ALLOCATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp_all=zeros(1,input.num_mes);
add_x_points=[];
bor_probe_spacing=1e9;
sur_probe_spacing=1e9;

tmp_A=[];
tmp_B=[];
tmp_M=[];
tmp_N=[];



pa=zeros(1,input.num_mes);
pb=zeros(1,input.num_mes);
pm=zeros(1,input.num_mes);
pn=zeros(1,input.num_mes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tmp_A(:,1)=input.ax;
tmp_A(:,2)=input.az;

tmp_B(:,1)=input.bx;
tmp_B(:,2)=input.bz;

tmp_M(:,1)=input.mx;
tmp_M(:,2)=input.mz;

tmp_N(:,1)=input.nx;
tmp_N(:,2)=input.nz;

% tmp_ALL matrix contains the UNIQUE x and y coordinates of all electrodes
% found
tmp_ALL=union(tmp_A,tmp_B,'rows');
tmp_ALL=union(tmp_ALL,tmp_M,'rows');
tmp_ALL=union(tmp_ALL,tmp_N,'rows');

clear tmp_A tmp_B tmp_M tmp_N


surface_electrodes=find(tmp_ALL(:,2)==0); % index of electrodes
borehole_electrodes=find(tmp_ALL(:,2)~=0); % index of electrodes


tmp_sur=tmp_ALL(surface_electrodes,:); % Electrodes with y coordinate =0
tmp_bor=tmp_ALL(borehole_electrodes,:); % All others



bor_pos=unique(tmp_bor(:,1)); % From all OTHERS, find different x-coordinates
mesh.num_of_bor=length(bor_pos); % number of boreholes 

sur_electrodes=unique(tmp_ALL(:,1)); % Keep all x-coordinates
bor_electrodes=unique(tmp_ALL(:,2)); % keep all y-coordinates

% Search if there are NO surface electrodes (compeare all x-coordiantes
% with number of boreholes
if length(sur_electrodes)==mesh.num_of_bor
    disp('No surface electrodes are found');
    surface_flag=1;
    sur_probe_spacing=1e9;
else
    disp('Surface elctrodes are found');
    surface_flag=0;
    % Find surface probe_spacing
    for i=1:length(sur_electrodes)-1
        tmp(i)=sur_electrodes(i+1,1)-sur_electrodes(i,1);
    end
    sur_probe_spacing=min(tmp); 
end


% Search for no boreholes 
if mesh.num_of_bor>0
    disp('Boreholes found');
    % Find borehole probe_spacing
    for i=1:length(bor_electrodes)-1
        tmp(i)=bor_electrodes(i+1,1)-bor_electrodes(i,1);
    end
    bor_probe_spacing=min(tmp); 
else
    disp('No boreholes found');
    bor_probe_spacing=1e9;
end


% find if borehoe on edge
if min(bor_pos) <= min(sur_electrodes)
    %flag to add some nodes on the edges
    add_bor_flag=1;
elseif max(bor_pos) >= max(sur_electrodes)
    %flag to add some nodes on the edges
    add_bor_flag=1;
else
    add_bor_flag=0;
end


% Now we may have br_min_probe_spacing and sr_probe_spacing. I only need
% one

mesh.probe_spacing=min(sur_probe_spacing,bor_probe_spacing);



% I must check again so add_x_pooint be excant between bor1 and bor2
if surface_flag==0
    add_x_points=sur_electrodes';
else 
    add_x_points=[min(bor_pos):mesh.probe_spacing:max(bor_pos)];
end


%ensure we have all borehole positions in add_x_points
add_x_points=[add_x_points bor_pos'];
add_x_points=unique(add_x_points);


% if we have boreholes on edge's add left and right more parameters
if add_bor_flag==1
    add_x_points=[min(add_x_points)-2*mesh.probe_spacing min(add_x_points)-mesh.probe_spacing ...
        add_x_points max(add_x_points)+mesh.probe_spacing max(add_x_points)+2*mesh.probe_spacing];
   add_x_points=unique(add_x_points);
end


% Create output matrix
mesh.add_x_points=add_x_points;

% Now make calculations comportable with Panos program
mesh.orig_probe_x=tmp_ALL(:,1);
mesh.orig_probe_z=tmp_ALL(:,2);
mesh.num_probes=length(mesh.orig_probe_x);


% Find probe coordiantes in numbering    
    for i=1:input.num_mes
      
        if input.m_array_type(i)==1      
            for j=1:mesh.num_probes
                if (input.ax(i)==mesh.orig_probe_x(j) && input.az(i)==mesh.orig_probe_z(j)) mesh.pa(i)=j; end
                if (input.bx(i)==mesh.orig_probe_x(j) && input.bz(i)==mesh.orig_probe_z(j)) mesh.pb(i)=j; end
                if (input.mx(i)==mesh.orig_probe_x(j) && input.mz(i)==mesh.orig_probe_z(j)) mesh.pm(i)=j; end
                if (input.nx(i)==mesh.orig_probe_x(j) && input.nz(i)==mesh.orig_probe_z(j)) mesh.pn(i)=j; end
            end
        end
      
        if input.m_array_type(i)==2      
            for j=1:mesh.num_probes
                if (input.ax(i)==mesh.orig_probe_x(j) && input.az(i)==mesh.orig_probe_z(j)) mesh.pa(i)=j;mesh.pb(i)=0; end                
                if (input.mx(i)==mesh.orig_probe_x(j) && input.mz(i)==mesh.orig_probe_z(j)) mesh.pm(i)=j; end
                if (input.nx(i)==mesh.orig_probe_x(j) && input.nz(i)==mesh.orig_probe_z(j)) mesh.pn(i)=j; end
            end
        end
     
        if input.m_array_type(i)==3      
            for j=1:mesh.num_probes
                if (input.ax(i)==mesh.orig_probe_x(j) && input.az(i)==mesh.orig_probe_z(j)) mesh.pa(i)=j;mesh.pb(i)=0; end
                if (input.mx(i)==mesh.orig_probe_x(j) && input.mz(i)==mesh.orig_probe_z(j)) mesh.pm(i)=j; mesh.pn(i)=0; end
            end
        end
 
    end

    
    
    % A check point to see if pa pb pm pn are all ok
    % Remember to check if zeros are found
    if length(pa)~=input.num_mes || length(pb)~=input.num_mes || length(pm)~=input.num_mes || length(pn)~=input.num_mes  
        disp('First error in measurment position');
        pause
    end
    
    for i=1:input.num_mes
        if input.m_array_type(i)==1
          tmp=min(mesh.pa);
          tmp=min(tmp,min(mesh.pb));
          tmp=min(tmp,min(mesh.pm));
          tmp=min(tmp,min(mesh.pn));
          if tmp==0
              disp('Something Went wrong');
              pause
          end
        elseif input.m_array_type(i)==3
          tmp=min(mesh.pa);
          tmp=min(tmp,min(mesh.pm));
          tmp=min(tmp,min(mesh.pn));
          if tmp==0
              disp('Something Went wrong');
              pause
          end
        elseif    input.m_array_type(i)==4
          tmp=min(mesh.pa);
          tmp=min(tmp,min(mesh.pm));
          if tmp==0
              disp('Something Went wrong');
              pause
          end
            
        end
        
    end    

end






function mesh=find_max_n(input,mesh)



mesh.depth_n=[];
% Find surface max_n base on A and M 
max_sr_n=max(sqrt((input.ax-input.mx).^2+(input.az-input.mz).^2));
max_sr_n=int32(max_sr_n/mesh.probe_spacing);



% If we have borehole, use electrode positions for layering
if  mesh.num_of_bor>0
    mesh.depth_n=unique(mesh.orig_probe_z);
    mesh.depth_n=[mesh.depth_n;0];
    mesh.depth_n=unique(mesh.depth_n);
    mesh.max_n=length(mesh.depth_n);
    tmp=1;
    for i=1:mesh.max_n-1
        new_depth_n(tmp)=mesh.depth_n(i);
        tmp=tmp+1;
        new_depth_n(tmp)=(mesh.depth_n(i)+mesh.depth_n(i+1))/2;
        tmp=tmp+1;
    end
    new_depth_n(tmp)=mesh.depth_n(mesh.max_n);
    
    %Here add one layer 
    new_depth_n(tmp+1)=new_depth_n(tmp)+mesh.probe_spacing;    
    
    mesh.depth_n=new_depth_n;
    mesh.max_n=length(mesh.depth_n);
    
else
    mesh.max_n=max_sr_n+1;
    if mesh.max_n>14
    mesh.max_n=14;
    end
%      mesh.max_n=9;
    max_br_n=0;
    mesh.depth_n=zeros(mesh.max_n,1);
    mesh.depth_n(1)=0;
    mesh.depth_n(2)=mesh.probe_spacing/2;
    m_factor=1.1; % This mean 10% thicker each layer
    for i=3:mesh.max_n
        tmp_thick=mesh.depth_n(i-1)-mesh.depth_n(i-2) ;
%         depth_n(i)=depth_n(i-1)+scale*probe_spacing/(2*probe_spacing);
          mesh.depth_n(i)=mesh.depth_n(i-1)+(m_factor)*tmp_thick;
    end
end

tmp=mesh.depth_n;
save('layers.mat','tmp');


mesh.depth_n=double(int32(mesh.depth_n*1000))/1000;

end

function mesh=custom_max_n(input,mesh)



% Custom Depth_n
% mesh.depth_n(1)=0;
% mesh.depth_n(2)=2.5;
% mesh.depth_n(3)=5.5;
% mesh.depth_n(4)=8.5;
% mesh.depth_n(5)=12;
% mesh.depth_n(6)=16;
% mesh.depth_n(7)=20;
% mesh.depth_n(8)=24;
% % mesh.depth_n(9)=29;
% % mesh.depth_n(10)=35;
% % mesh.depth_n(11)=41;
% % mesh.depth_n(12)=15;
% % mesh.depth_n(13)=15.758;
% % mesh.depth_n(14)=18.234;
% % mesh.depth_n(15)=20.958;
% % % 
% 
% mesh.max_n=8;
% mesh.depth_n=mesh.depth_n(1:mesh.max_n);
% 


mesh.depth_n=double(int32(mesh.depth_n*1000))/1000;

end

function [mesh]=calc_param(input,mesh)



L=zeros(mesh.max_n,1);
if mesh.num_of_bor==0
    L(1)=0;
    for i=2:mesh.max_n
        j=int32(mod(double(i),2.0));
        
        if i<4 
            j=1;
        elseif i>=4 && i<7
            j=2;
        elseif i>=7 
            j=3;
        end
        
        j=1;
%         j=0;
        L(i)=L(i-1)+j;
    end
end

if mesh.num_of_bor==2
    L=zeros(mesh.max_n,1);
end
 L=zeros(mesh.max_n,1);

scale=35.5*mesh.probe_spacing/4;
left_boundary=min(mesh.add_x_points)-mesh.probe_spacing*scale;
right_boundary=max(mesh.add_x_points)+mesh.probe_spacing*scale;
bottom=max(mesh.depth_n)+mesh.probe_spacing*scale;

tmp1=0;
for i=1:mesh.max_n-1
       
        map=0;
                
            d1=mesh.depth_n(i);
            d2=mesh.depth_n(i+1);
               
        for J=L(i)+1:length(mesh.add_x_points)-L(i)-1
            
          tmp1=tmp1+1;
		  map=map+1;
		  mesh.map_param(i,map)=tmp1;
          

      mesh.tmp_param(tmp1,3)=mesh.add_x_points(J);  %x1
	  mesh.tmp_param(tmp1,4)=mesh.add_x_points(J+1);  %x2
                    
      mesh.tmp_param(tmp1,5)=d1;                         %y1
      mesh.tmp_param(tmp1,6)=d2;                         %y2


      if(J==L(i)+1) ;mesh.tmp_param(tmp1,3)=min(mesh.add_x_points)-mesh.probe_spacing;end                           % left
	  if(J==length(mesh.add_x_points)-L(i)-1) ;mesh.tmp_param(tmp1,4)=max(mesh.add_x_points)+mesh.probe_spacing;end % right 
          
      mesh.tmp_param(tmp1,1)=(mesh.tmp_param(tmp1,3)+mesh.tmp_param(tmp1,4)) /2; % x center
	  mesh.tmp_param(tmp1,2)=(d1+d2)/2;                 % y center

      
      
      
      
	   mesh.param_x(tmp1,1)=(mesh.add_x_points(J)+mesh.add_x_points(J+1)) /2;
	   mesh.param_y(tmp1,1)=(d1+d2)/2;

   %    if(dp>depth_n_max) depth_n_max=dp;end
         end

            
        end
    
    
  mesh.num_param=tmp1;  

  %clear x_scaled
end


function mesh=create_mesh_gen(input,mesh)






mesh.no_elm_per_param=zeros(mesh.num_param,1);


node=[];
edge=[];
faces=[];
icon=[];
x_node_coord=[];
y_node_coord=[];
mes_nodes=[];


tmp=1;
for i=1:mesh.num_param
    
    node(tmp,1)=mesh.tmp_param(i,3); % x1
    node(tmp,2)=mesh.tmp_param(i,5); % y1
    tmp=tmp+1;
    node(tmp,1)=mesh.tmp_param(i,3); % x1
    node(tmp,2)=mesh.tmp_param(i,6); % y2
    tmp=tmp+1;
    node(tmp,1)=mesh.tmp_param(i,4); % x2
    node(tmp,2)=mesh.tmp_param(i,5); % y1
    tmp=tmp+1;
    node(tmp,1)=mesh.tmp_param(i,4); % x2
    node(tmp,2)=mesh.tmp_param(i,6); % y2
    tmp=tmp+1;
    
    % Interpolation points (i need the to add more elements within each
    % parameter
    node(tmp,1)=(mesh.tmp_param(i,3)+mesh.tmp_param(i,4))/2; %x1y1---x2y1
    node(tmp,2)=mesh.tmp_param(i,5);
    tmp=tmp+1;
    
    node(tmp,1)=mesh.tmp_param(i,4); %x2y1---x2y2
    node(tmp,2)=(mesh.tmp_param(i,5)+mesh.tmp_param(i,6))/2;
    tmp=tmp+1;
    
    node(tmp,1)=(mesh.tmp_param(i,3)+mesh.tmp_param(i,4))/2; %x2y2---x1y2
    node(tmp,2)=mesh.tmp_param(i,6);
    tmp=tmp+1;    
    
    node(tmp,1)=mesh.tmp_param(i,3); %x1y2---x1y1
    node(tmp,2)=(mesh.tmp_param(i,5)+mesh.tmp_param(i,6))/2;
    tmp=tmp+1;    
    
   
end


% Make sure we have nodes on the electrode positions

tmp_node(:,1)=mesh.orig_probe_x;
tmp_node(:,2)=mesh.orig_probe_z;
% remove duplicate points
node=union(node,tmp_node,'rows'); 

clear tmp_node
%Here add node_coordinates for boundary conditions
%    (1)--------------(2)                    (3)-----------(4)
%     |                |                      |             |
%     |                |                      |             |
%     |                |                      |             |
%     |               (5)--------------------(6)            |
%     |                                                     |
%     |                                                     |
%    (7)---------------------------------------------------(8)


scale=35.5*mesh.probe_spacing/2;
left_boundary=min(mesh.add_x_points)-mesh.probe_spacing*scale;
right_boundary=max(mesh.add_x_points)+mesh.probe_spacing*scale;
bottom=max(mesh.depth_n)+mesh.probe_spacing*scale;


%1
tmp_node(1,1)=left_boundary; 
tmp_node(1,2)=0;

%2
tmp_node(2,1)=min(mesh.tmp_param(:,3)); 
tmp_node(2,2)=0;

%3
tmp_node(3,1)=max(mesh.tmp_param(:,4)); 
tmp_node(3,2)=0;

%4
tmp_node(4,1)=right_boundary; 
tmp_node(4,2)=0;

%5
tmp_node(5,1)=min(mesh.tmp_param(:,3)); 
tmp_node(5,2)=max(mesh.tmp_param(:,6));

%6
tmp_node(6,1)=max(mesh.tmp_param(:,4));  
tmp_node(6,2)=max(mesh.tmp_param(:,6));

%7
tmp_node(7,1)=left_boundary;  
tmp_node(7,2)=bottom;

%8
tmp_node(8,1)=right_boundary;  
tmp_node(8,2)=bottom;

% remove duplicate points
node=union(node,tmp_node,'rows'); 






tmp=1;
tmp2=1;
%now create the edges

edge=zeros(4*mesh.num_param,2);
for i=1:mesh.num_param+1
    faces{i}=[];
end



for i=1:mesh.num_param
    
    x1=mesh.tmp_param(i,3);
    x2=mesh.tmp_param(i,4);
    y1=mesh.tmp_param(i,5);
    y2=mesh.tmp_param(i,6);
    
    % x1,y1 --> x2,y1
            
            ind_all_x=find(node(:,1)<=x2 & node(:,2)==y1);
            ind_all_x2=find(node(:,1)>=x1 & node(:,2)==y1);
            ind_all_x=intersect(ind_all_x,ind_all_x2);
            
            
            for j=1:length(ind_all_x)-1
                ind_x1=ind_all_x(j);
                ind_x2=ind_all_x(j+1);
                edge(tmp,1)=ind_x1;
                edge(tmp,2)=ind_x2;
                tmp=tmp+1;  
                faces{i}=[faces{i} tmp2];
                tmp2=tmp2+1;
            end
        
        
    % x2,y1 --> x2,y2
            
            ind_all_x=find(node(:,1)==x2 & node(:,2)>=y1);
            ind_all_x2=find(node(:,1)==x2 & node(:,2)<=y2);
            ind_all_x=intersect(ind_all_x,ind_all_x2);
            
            
            for j=1:length(ind_all_x)-1
                ind_x1=ind_all_x(j);
                ind_x2=ind_all_x(j+1);
                edge(tmp,1)=ind_x1;
                edge(tmp,2)=ind_x2;
                tmp=tmp+1;  
                faces{i}=[faces{i} tmp2];
                tmp2=tmp2+1;
            end
        
     % x2,y2 --> x1,y2
            
            ind_all_x=find(node(:,1)>=x1 & node(:,2)==y2);
            ind_all_x2=find(node(:,1)<=x2 & node(:,2)==y2);
            ind_all_x=intersect(ind_all_x,ind_all_x2);
            
            
            for j=1:length(ind_all_x)-1
                ind_x1=ind_all_x(j);
                ind_x2=ind_all_x(j+1);
                edge(tmp,1)=ind_x1;
                edge(tmp,2)=ind_x2;
                tmp=tmp+1;  
                faces{i}=[faces{i} tmp2];
                tmp2=tmp2+1;
            end
        
     % x1,y2 --> x1,y1
            
            ind_all_x=find(node(:,1)==x1 & node(:,2)>=y1);
            ind_all_x2=find(node(:,1)==x1 & node(:,2)<=y2);
            ind_all_x=intersect(ind_all_x,ind_all_x2);
            
            
            for j=1:length(ind_all_x)-1
                ind_x1=ind_all_x(j);
                ind_x2=ind_all_x(j+1);
                edge(tmp,1)=ind_x1;
                edge(tmp,2)=ind_x2;
                tmp=tmp+1;  
                faces{i}=[faces{i} tmp2];
                tmp2=tmp2+1;
            end        
 
end


% Here create face for the boundary conditions
%1->2
ind_x1=find(node(:,1)==left_boundary & node(:,2)==0);
ind_x2=find(node(:,1)==min(mesh.tmp_param(:,3)) & node(:,2)==0);
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;

%2->5
ind_x1=find(node(:,1)==min(mesh.tmp_param(:,3)) & node(:,2)==0);
ind_x2=find(node(:,1)==min(mesh.tmp_param(:,3)) & node(:,2)==max(mesh.tmp_param(:,6)));
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;


%5->6
ind_x1=find(node(:,1)==min(mesh.tmp_param(:,3)) & node(:,2)==max(mesh.tmp_param(:,6)));
ind_x2=find(node(:,1)==max(mesh.tmp_param(:,4)) & node(:,2)==max(mesh.tmp_param(:,6)));
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;


%6->3
ind_x1=find(node(:,1)==max(mesh.tmp_param(:,4)) & node(:,2)==max(mesh.tmp_param(:,6)));
ind_x2=find(node(:,1)==max(mesh.tmp_param(:,4)) & node(:,2)==0);
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;



%3->4
ind_x1=find(node(:,1)==max(mesh.tmp_param(:,4)) & node(:,2)==0);
ind_x2=find(node(:,1)==right_boundary & node(:,2)==0);
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;



%4->8
ind_x1=find(node(:,1)==right_boundary & node(:,2)==0);
ind_x2=find(node(:,1)==right_boundary & node(:,2)==bottom);
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;



%8->7
ind_x1=find(node(:,1)==right_boundary & node(:,2)==bottom);
ind_x2=find(node(:,1)==left_boundary & node(:,2)==bottom);
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;



%7->1
ind_x1=find(node(:,1)==left_boundary & node(:,2)==bottom);
ind_x2=find(node(:,1)==left_boundary & node(:,2)==0);
edge(tmp,1)=ind_x1;
edge(tmp,2)=ind_x2;
tmp=tmp+1;
faces{i+1}=[faces{i+1} tmp2];
tmp2=tmp2+1;



orig_probe_x=unique(node(:,1));
save('orig_probe_x.mat','orig_probe_x');
anaglyfo;

h=gcf;
uiwait(h);
load('anaglyfo_data.mat');
% % % % topography here

% % % % Include topography in orig_probe_x
% % % 
for i=1:length(anaglyfo_data(:,1))
    for j=1:length(mesh.orig_probe_x)
        if mesh.orig_probe_x(j)==anaglyfo_data(i,1)
            mesh.orig_probe_z(j)=mesh.orig_probe_z(j)-anaglyfo_data(i,2);
        end
    end
end
tmp=unique(node(:,2));

for i=1:length(node(:,1))
    for j=1:length(anaglyfo_data)
         if node(i,1)==anaglyfo_data(j,1) && node(i,2)~=tmp(end-1)
            node(i,2)=node(i,2)-anaglyfo_data(j,2);
        end
    end
end






clear gca


plot(mesh.orig_probe_x,-mesh.orig_probe_z,'o');
hold on
%Here I will create specific element stize on each edge. I am searching for
%edges that are on theleft, right and bottom, and all the rest



h = waitbar(0,'Mesh Design.Please wait up to minutes....');
hdata.hmax = mesh.probe_spacing/3;
% hdata.edgeh=ll;
options.output=false;
warning off
addpath ('mesh2d');
[p ,t,fnum]=meshfaces(node,edge,faces,[],options);
% [p ,t,fnum]=meshfaces(node,edge,faces,hdata,options);
waitbar(1/1);
close(h);
warning on


mesh.faces=faces;
mesh.node=node;
mesh.edge=edge;

% And also correct the x,y centers of parameters
% % %  for i=1:mesh.num_param
% % %     tmp_faces=faces{i};
% % %     tmp_edges=edge(tmp_faces);
% % %     %x1=min(node(tmp_edges,1))*probe_spacing/2;
% % %     %x2=max(node(tmp_edges,1))*probe_spacing/2;
% % %     y1=min((node(tmp_edges,2)))*mesh.probe_spacing/2;
% % %     y2=max((node(tmp_edges,2)))*mesh.probe_spacing/2;
% % %     
% % %     %param_x(i)=(x1+x2)/2;
% % %     %param_y(i)=(y1+y2)/2;    
% % %  end


mesh.num_nodes=length(p(:,2));
mesh.num_elements=length(fnum);

mesh.x_node_coord=p(:,1);
mesh.y_node_coord=p(:,2);
mesh.icon(1,:)=t(:,1);
mesh.icon(2,:)=t(:,2);
mesh.icon(3,:)=t(:,3);
mesh.icon(4,:)=fnum(:);









p(:,2)=-p(:,2);
node(:,2)=-node(:,2);

%figure('Name','Mesh')
  % plot(x_node_coord,y_node_coord,'b.','markersize',1)
   hold on;
   % Colour mesh for each face
   col = ['b','r','g','k','m'];
   for k = 1:length(faces)
      colk = mod(k,length(col));
      if (colk==0)
         colk = length(col);
      end
       patch('faces',t(fnum==k,:),'vertices',p,'facecolor','w','edgecolor',col(colk));
   end
   patch('faces',edge,'vertices',node,'facecolor','none','edgecolor','k')
   % Highlisght low q triangles in debug mode
 %  if options.debug
 %     pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3.0;
 %     plot(pc(q<0.5,1),pc(q<0.5,2),'r.')
 %  end
   axis equal;% off;







%/* Find probe coordinates */

for i=1:mesh.num_probes
    for j=1:mesh.num_nodes
        if( (mesh.x_node_coord(j)==mesh.orig_probe_x(i)) && (mesh.y_node_coord(j)==mesh.orig_probe_z(i)))
            mesh.mes_nodes(i)=j;
        end
    end
end
% Find a electrode for sink


    for j=1:mesh.num_nodes
        if( (mesh.x_node_coord(j)<0.8*min(mesh.orig_probe_x) && (mesh.y_node_coord(j)==0)))
            mesh.sink=j;
        end
    end



% Second check 
if length(mesh.mes_nodes)~=mesh.num_probes || min(mesh.mes_nodes)==0
    disp('Second error in measurements nodes');
    pause
end


%/* ******* FIND EBC NODES ********* */

left_ebc=find(mesh.x_node_coord(:)==min(mesh.x_node_coord));
right_ebc=find(mesh.x_node_coord(:)==max(mesh.x_node_coord));
bottom_ebc=find(mesh.y_node_coord(:)==max(mesh.y_node_coord));

mesh.node_ebc=union(left_ebc,right_ebc);
mesh.node_ebc=union(mesh.node_ebc,bottom_ebc);

mesh.num_ebc=length(mesh.node_ebc);



for i=1:mesh.num_param
    k=0;
    for j=1:mesh.num_elements
        if(mesh.icon(4,j)==i)
            k=k+1;
            mesh.param(i,k)=j;
        end
    end
end


% Find number of elements in each parameter
for j=1:mesh.num_param
    [la la2]=find(mesh.param(j,:)~=0);
    mesh.no_elm_per_param(j)=length(la2);
end
    
  

hold off
end


	

%--------------------------------------------------------------------------
%                         MESH GRID DATA                                  |
%--------------------------------------------------------------------------

function mesh=create_mesh_grid(input,mesh)

% Add topography
load('anaglyfo_data.mat');
y_tmp=zeros(mesh.num_param,1);
for i=1:mesh.num_param
    %check for 1st parameter
    if mesh.param_x(i)==min(mesh.param_x)
        x_tmp=anaglyfo_data(3,1);
        y_tmp(i)=anaglyfo_data(3,2);
    elseif mesh.param_x(i)==max(mesh.param_x)
        x_tmp=anaglyfo_data(end-2,1);
        y_tmp(i)=anaglyfo_data(3,2);
    else
        x_tmp=mesh.param_x(i);
        ind=find(anaglyfo_data(:,1)==x_tmp);
        y_tmp(i)=anaglyfo_data(ind,2);
    end
    
    
    
end


mesh.min_x=min(mesh.param_x);
mesh.max_x=max(mesh.param_x);

mesh.min_y=min(mesh.param_y-y_tmp);
mesh.max_y=max(mesh.param_y-y_tmp);


xx=union(mesh.param_x, mesh.param_x); % take the unique x_line
yy=union(mesh.param_y-y_tmp, mesh.param_y-y_tmp); %take the unique y_line

tmp_x=xx(2)-xx(1);% find the x step
tmp_y=yy(2)-yy(1);% find the y_step
yy=-yy; %Now creative y as negative



[mesh.XI mesh.YI]=meshgrid(xx,yy);% and now create the gridata 



xnew=[mesh.min_x:tmp_x/2:mesh.max_x];
ynew=[mesh.min_y:tmp_y/2:mesh.max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[mesh.xxx,mesh.yyy]=meshgrid(xnew,ynew);



mesh.y_tmp=y_tmp;
% if bgr_res_flag==2
%     
%     figure;
%     ZI=griddata(param_x-neg_x,-param_y,log10(bgr_param),XI,YI);
%     
%     ZII = interp2(XI,YI,ZI,xxx,yyy,'cubic');
%     pcolor(xxx,yyy,ZII);
%     shading interp
%     axis equal
%     title('Forward model');
%     xlim([min_x max_x]);
%     ylim([-max_y -min_y]);
%     xlabel('Distance (m)');
%     ylabel('Depth (m)');
%     colorbar;
%     drawnow;
%     
% end



end

