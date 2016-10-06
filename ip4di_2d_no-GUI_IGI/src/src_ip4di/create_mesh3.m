%----------------------------------------------------------------------%
%
% Function CREATE_MESH3: creates the finite-element mesh.
%
% Author: Marios Karaoulis, Colorado School of Mines
%
%----------------------------------------------------------------------%

% Modified by Francois Lavoue, Colorado School of Mines, September 2015.
%   Modifications include
%   1) some comments for more readability.
%   2) the possibility to load a pre-existing mesh.
%   3) the externalization of plot routines (see plot_mesh.m).
%   4) the externalization of the definition of mean_res (see define_mean_res.m),
%      because mean_res must match the model background value for the solution to be accurate.
%   5) refine the mesh between electrodes and extend it laterally.
%   6) deal with some numerical issues in the removal of duplicate points
%      (use the manual criterion input.accuracy_threshold instead of == boolean).
%   7) avoid the call to the anaglyfo function, to suppress the need for user's manual intervention.

% Modified by Francois Lavoue, August 2016.
%   -> re-introduce topography, following the original approach of anaglyfo function,
%      ie 1- electrode coordinates given in input data file must NOT contain topography,
%            such that surface and borehole electrodes can be determined.
%         2- then topography is read from file, interpolated at electrode and node locations,
%            and used to correct electrode and node z-coordinates.
%   -> change variable name 'mesh.anaglyfo_data' into 'mesh.topo_nodes'.
%   -> introduce new vector mesh.topo_param to store interpolated topographic correction
%      (useful for further plots).

function [mesh]=create_mesh3(input,mesh)

disp(' ')
disp('-------------------------')
disp('    ENTER CREATE_MESH    ')
disp('-------------------------')
disp(' ')

if input.read_mesh==1
   % READ MESH
   mesh=load(input.file_mesh_in);
   disp(sprintf('Read previously generated mesh in file %s', input.file_mesh_in))

else   % BUILD MESH

   t1=cputime;

   mesh=find_probe_spacing(input);
   if input.depth_n==0;
      mesh=find_max_n(input,mesh);
   else
      mesh.max_n=length(input.depth_n);
      mesh.depth_n=input.depth_n;
   end
   mesh=calc_param(input,mesh);
   mesh=create_mesh_gen(input,mesh);

   % calculate_half_band
   mesh=create_mesh_grid(input,mesh);
   mesh=geometrical_factor(input,mesh);

   t2=cputime;
   disp(sprintf('TIME FOR GENERATING MESH = %f',t2-t1));

end   %end if re-use mesh


% save mesh
try
    save(input.file_mesh_out, '-struct', 'mesh');
catch
    msgbox('cannot save mesh');
end   %end try


% plot mesh
if input.plot_mesh==1
   plot_mesh(mesh);
end


%% define mean res. (FL: now done in define_mean_res.m called in main program)
%mesh=define_mean_res(input,mesh);

%% FL: lines below are now included in define_mean_res.m
%% /!\ Default values must match homogeneous background,
%%     otherwise solution is wrong
%mesh.init_res_param(1:mesh.num_param)=mesh.mean_res;
%mesh.res_param1(1:mesh.num_param,1)=mesh.mean_res;
%mesh.res_param2(1:mesh.num_param,1)=mesh.mean_res;
%mesh.prop(1:mesh.num_elements)=mesh.mean_res;
%
%if input.time_lapse_flag==1
%    mesh.d4_res_param1=mesh.mean_res*ones(mesh.num_param,input.num_files);
%    mesh.d4_res_param2=mesh.mean_res*ones(mesh.num_param,input.num_files);   
%end

end   %-- end function create_mesh
%----------------------------------------------------------------------%



%----------------------------------------------------------------------%
%
%   Function FIND_PROBE_SPACING: find electrode locations and define
%                                mesh nodes accordingly.
%
% The main idea for this function is to calculate the basic dimensions
% of the mesh. I need to find the surface electrodes, the borehole elec-
% trodes, the positions of the electrodes etc. If for some reason a bo-
% rehole is on edge or if there are no surface electrodes I must ensure
% to add the neccessary points.
%
%----------------------------------------------------------------------%

function mesh=find_probe_spacing(input)

%-- pre-allocation
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
%--

tmp_A(:,1)=input.ax;
tmp_A(:,2)=input.az;

tmp_B(:,1)=input.bx;
tmp_B(:,2)=input.bz;

tmp_M(:,1)=input.mx;
tmp_M(:,2)=input.mz;

tmp_N(:,1)=input.nx;
tmp_N(:,2)=input.nz;

% tmp_ALL matrix contains the UNIQUE x and y coordinates of all electrodes found
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


% Search if there are NO surface electrodes
% (compare all x-coordinates with number of boreholes)
if length(sur_electrodes)==mesh.num_of_bor
    disp('No surface electrodes are found');
    surface_flag=1;
    sur_probe_spacing=1e9;
else
    disp(sprintf('%i Surface electrodes are found',length(sur_electrodes)));
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


% Find if borehole must be added on edges
if min(bor_pos) <= min(sur_electrodes)
   % flag to add some nodes on the edges
   add_bor_flag=1;
elseif max(bor_pos) >= max(sur_electrodes)
   %flag to add some nodes on the edges
   add_bor_flag=1;
else
   add_bor_flag=0;
end


% Now we may have bor_probe_spacing and sur_probe_spacing. I only need one
mesh.probe_spacing=min(sur_probe_spacing,bor_probe_spacing);


% I must check again so add_x_point be excant between bor1 and bor2
if surface_flag==0
% in case of surface electrodes, refine and extend laterally
    add_x_points=[min(sur_electrodes)-input.x_margin : mesh.probe_spacing/input.ncells_per_elec : max(sur_electrodes)+input.x_margin];
    %add_x_points=sur_electrodes';   %FL commented
else   % no surface electrode: do not extend laterally, this is done below
    add_x_points=[min(bor_pos):mesh.probe_spacing/input.ncells_per_elec:max(bor_pos)];
end

% ensure we have all borehole positions in add_x_points
add_x_points=[add_x_points bor_pos'];

% round before removing duplicate points,
% otherwise duplicate points may be considered as different due to rounding errors
add_x_points=round(add_x_points/input.accuracy_threshold)*input.accuracy_threshold;

% remove duplicate points
add_x_points=unique(add_x_points);

if add_bor_flag==1
% if we have boreholes on edges, add two parameter cells on left and right
   add_x_points=[ min(add_x_points)-2*mesh.probe_spacing min(add_x_points)-mesh.probe_spacing ...
                  add_x_points ...
                  max(add_x_points)+mesh.probe_spacing max(add_x_points)+2*mesh.probe_spacing ];
   add_x_points=unique(add_x_points);
else   % add just one parameter cell on left and right
   add_x_points=[ min(add_x_points)-mesh.probe_spacing ...
                  add_x_points ...
                  max(add_x_points)+mesh.probe_spacing ];
   add_x_points=unique(add_x_points);
end

% Create output matrix
mesh.add_x_points=add_x_points;

% Now make calculations comportable with Panos program
mesh.orig_probe_x=tmp_ALL(:,1);
mesh.orig_probe_z=tmp_ALL(:,2);
mesh.num_probes=length(mesh.orig_probe_x);


 %-- Find probe coordinates in numbering
 for i=1:input.num_mes

     if input.m_array_type(i)==1      
        for j=1:mesh.num_probes
            if (input.ax(i)==mesh.orig_probe_x(j) && input.az(i)==mesh.orig_probe_z(j)) mesh.pa(i)=j; end
            if (input.bx(i)==mesh.orig_probe_x(j) && input.bz(i)==mesh.orig_probe_z(j)) mesh.pb(i)=j; end
            if (input.mx(i)==mesh.orig_probe_x(j) && input.mz(i)==mesh.orig_probe_z(j)) mesh.pm(i)=j; end
            if (input.nx(i)==mesh.orig_probe_x(j) && input.nz(i)==mesh.orig_probe_z(j)) mesh.pn(i)=j; end

            %if ( abs(input.ax(i)-mesh.orig_probe_x(j))<input.accuracy_threshold && ...
            %     abs(input.az(i)-mesh.orig_probe_z(j))<input.accuracy_threshold )  ...
            %   mesh.pa(i)=j; end
            %if ( abs(input.bx(i)-mesh.orig_probe_x(j))<input.accuracy_threshold && ...
            %     abs(input.bz(i)-mesh.orig_probe_z(j))<input.accuracy_threshold )  ...
            %   mesh.pb(i)=j; end
            %if ( abs(input.mx(i)-mesh.orig_probe_x(j))<input.accuracy_threshold && ...
            %     abs(input.mz(i)-mesh.orig_probe_z(j))<input.accuracy_threshold )  ...
            %   mesh.pm(i)=j; end

            %if ( abs(input.nx(i)-mesh.orig_probe_x(j))<input.accuracy_threshold && ...
            %     abs(input.nz(i)-mesh.orig_probe_z(j))<input.accuracy_threshold )  ...
            %   mesh.pn(i)=j; end
        end
     end

     if input.m_array_type(i)==2      
        for j=1:mesh.num_probes
            if (input.ax(i)==mesh.orig_probe_x(j) && input.az(i)==mesh.orig_probe_z(j)) mesh.pa(i)=j; mesh.pb(i)=0; end                
            if (input.mx(i)==mesh.orig_probe_x(j) && input.mz(i)==mesh.orig_probe_z(j)) mesh.pm(i)=j; end
            if (input.nx(i)==mesh.orig_probe_x(j) && input.nz(i)==mesh.orig_probe_z(j)) mesh.pn(i)=j; end
        end
     end

     if input.m_array_type(i)==3      
        for j=1:mesh.num_probes
            if (input.ax(i)==mesh.orig_probe_x(j) && input.az(i)==mesh.orig_probe_z(j)) mesh.pa(i)=j; mesh.pb(i)=0; end
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

 end   %end for i=1:input.num_mes

end   %end function mesh=find_probe_spacing(input)
%----------------------------------------------------------------------%



%----------------------------------------------------------------------%
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

if input.debug==1
   % output layers
   tmp=mesh.depth_n;
   save('layers.mat','tmp');
end

mesh.depth_n=double(int32(mesh.depth_n*1000))/1000;

end   %end function mesh=find_max_n(input,mesh)
%----------------------------------------------------------------------%



%----------------------------------------------------------------------%
%%FL: commented because not used
%function mesh=custom_max_n(input,mesh)
%
%% Custom Depth_n
%% mesh.depth_n(1)=0;
%% mesh.depth_n(2)=2.5;
%% mesh.depth_n(3)=5.5;
%% mesh.depth_n(4)=8.5;
%% mesh.depth_n(5)=12;
%% mesh.depth_n(6)=16;
%% mesh.depth_n(7)=20;
%% mesh.depth_n(8)=24;
%% % mesh.depth_n(9)=29;
%% % mesh.depth_n(10)=35;
%% % mesh.depth_n(11)=41;
%% % mesh.depth_n(12)=15;
%% % mesh.depth_n(13)=15.758;
%% % mesh.depth_n(14)=18.234;
%% % mesh.depth_n(15)=20.958;
%% % % 
%% 
%% mesh.max_n=8;
%% mesh.depth_n=mesh.depth_n(1:mesh.max_n);
%% 
%
%mesh.depth_n=double(int32(mesh.depth_n*1000))/1000;
%
%end
%----------------------------------------------------------------------%



%----------------------------------------------------------------------%
%
%   Function CALC_PARAM: define grid of parameters cells based on elec-
%                        trode locations.
%
%----------------------------------------------------------------------%

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


%- Define boundaries
if input.user_boundary_flag==1
   % define boundaries according to user-defined PML width
   left_boundary=min(mesh.add_x_points)-input.extra_width;
   right_boundary=max(mesh.add_x_points)+input.extra_width;
   bottom=max(mesh.depth_n)+input.extra_width;

elseif input.user_boundary_flag==2
   % variation of default boundaries using Marios' scale
   scale=35.5*mesh.probe_spacing/4;
   left_boundary=min(mesh.add_x_points)-scale;   % do NOT reuse probe_spacing here,
   right_boundary=max(mesh.add_x_points)+scale;  % such that the added value has well the dimension of spacing
   bottom=max(mesh.depth_n)+scale;

else
   % default boundaries using Marios' scale
   scale=35.5*mesh.probe_spacing/4;    % /!\ HERE IT IS 4, BELOW IT IS 2???
   %scale=35.5*mesh.probe_spacing/2;
   left_boundary=min(mesh.add_x_points)-mesh.probe_spacing*scale;
   right_boundary=max(mesh.add_x_points)+mesh.probe_spacing*scale;
   bottom=max(mesh.depth_n)+mesh.probe_spacing*scale;
end


%-- Define coordinates of parameter cells
tmp1=0;
% loop over layers
for i=1:mesh.max_n-1

    map=0;

    d1=mesh.depth_n(i);
    d2=mesh.depth_n(i+1);

    for J=L(i)+1:length(mesh.add_x_points)-L(i)-1

        tmp1=tmp1+1;
	map=map+1;
	mesh.map_param(i,map)=tmp1;

        mesh.tmp_param(tmp1,3)=mesh.add_x_points(J);    %x1
        mesh.tmp_param(tmp1,4)=mesh.add_x_points(J+1);  %x2

        mesh.tmp_param(tmp1,5)=d1;                      %y1
        mesh.tmp_param(tmp1,6)=d2;                      %y2

        %-- put a margin of width=probe_spacing to the left and to the right of the model
        %FL: Commented because it modifies the cell defined above and create a cell that has twice the width of other cells.
        %    This can cause problems in case of topography: as the mesher will link the vertices of this extra-width cell,
        %    the 1st and last electrodes located between the vertices may not be on the surface anymore.
        %    Instead, additional x-coords have been added to mesh.add_x_points in find_probe_spacing,
        %    and they are treated just like other x-points here, resulting in additional parameter cells on left and right.
        %if(J==L(i)+1); mesh.tmp_param(tmp1,3)=min(mesh.add_x_points)-mesh.probe_spacing; end                           % left
        %if(J==length(mesh.add_x_points)-L(i)-1); mesh.tmp_param(tmp1,4)=max(mesh.add_x_points)+mesh.probe_spacing; end % right 

        mesh.tmp_param(tmp1,1)=(mesh.tmp_param(tmp1,3)+mesh.tmp_param(tmp1,4)) /2; % x center
        mesh.tmp_param(tmp1,2)=(d1+d2)/2;                                          % y center

        mesh.param_x(tmp1,1)=(mesh.add_x_points(J)+mesh.add_x_points(J+1)) /2;
        mesh.param_y(tmp1,1)=(d1+d2)/2;

        %if(dp>depth_n_max) depth_n_max=dp;end
     end   %end     for J=L(i)+1:length(mesh.add_x_points)-L(i)-1

end   %-- end loop over layers (for i=1:mesh.max_n-1)

%-- define nb of parameter (cells)
mesh.num_param=tmp1;  

end   %-- end function calc_param
%----------------------------------------------------------------------%



%----------------------------------------------------------------------%
%
%   Function CREATE_MESH_GEN: generate finite-element mesh based on
%                             parameter grid.
%
%----------------------------------------------------------------------%

function mesh=create_mesh_gen(input,mesh)

% init.
tmp=1;
node=[];
edge=[];
faces=[];
icon=[];
x_node_coord=[];
y_node_coord=[];
mesh.mes_nodes=[];
mesh.no_elm_per_param=zeros(mesh.num_param,1);


%-- loop over parameter cells
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

    % Interpolation points (add more elements within each parameter cell)
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


%-- Make sure we have nodes on the electrode positions
tmp_node(:,1)=mesh.orig_probe_x;
tmp_node(:,2)=mesh.orig_probe_z;
% remove duplicate points
node=union(node,tmp_node,'rows'); 

clear tmp_node


%-- Add node_coordinates for boundary conditions
%
%    (1)--------------(2)                    (3)-----------(4)
%     |                |                      |             |
%     |                |                      |             |
%     |                |                      |             |
%     |               (5)--------------------(6)            |
%     |                                                     |
%     |                                                     |
%    (7)---------------------------------------------------(8)

%- define boundaries
if input.user_boundary_flag==1
   % define boundaries according to user-defined PML width
   left_boundary=min(mesh.add_x_points)-input.extra_width;
   right_boundary=max(mesh.add_x_points)+input.extra_width;
   bottom=max(mesh.depth_n)+input.extra_width;

elseif input.user_boundary_flag==2
   % variation of default boundaries using Marios' scale
   scale=35.5*mesh.probe_spacing/2;
   left_boundary=min(mesh.add_x_points)-scale;   % do NOT reuse probe_spacing here,
   right_boundary=max(mesh.add_x_points)+scale;  % such that the added value has well the dimension of spacing
   bottom=max(mesh.depth_n)+scale;

else
   % default boundaries using Marios' scale
   scale=35.5*mesh.probe_spacing/2;    %FL /!\ HERE IT IS 2, BEFORE IT WAS 4???
   left_boundary=min(mesh.add_x_points)-mesh.probe_spacing*scale;
   right_boundary=max(mesh.add_x_points)+mesh.probe_spacing*scale;
   bottom=max(mesh.depth_n)+mesh.probe_spacing*scale;
end


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

% remove similar points up to accuracy threshold manually
iduplicate=[];
for inode=1:length(node)
    ifind=find( abs(node(:,1)-node(inode,1))<input.accuracy_threshold ...
              & abs(node(:,2)-node(inode,2))<input.accuracy_threshold );
    if length(ifind)>1
       iduplicate=[iduplicate ; ifind(2:end)];
    end
end
node(iduplicate,:)=[];



%-- Create edges of parameter cells

% init.
tmp=1;
tmp2=1;
edge=zeros(4*mesh.num_param,2);
for i=1:mesh.num_param+1
    faces{i}=[];
end

%-- loop over parameter cells
for i=1:mesh.num_param

    % cell vertices coordinates
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

end   %-- end loop over parameter cells (for i=1:mesh.num_param)


%-- Create faces for the boundary conditions

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


%-- Topography
% init: no topo
topo_nodes(:,1)=unique(node(:,1));
topo_nodes(:,2)=0;

%-- add topo
if input.topo_flag==1
   %- read topo file
   topo=load(input.file_topo_in);

   %- add topo values if missing at model ends
   if topo_nodes(1,1) < topo(1,1)
      topo=[ [topo_nodes(1,1);topo(:,1)] [topo(1,2);topo(:,2)] ];
   end
   if topo_nodes(end,1) > topo(end,1)
      topo=[ [topo(:,1);topo_nodes(end,1)] [topo(:,2);topo(end,2)] ];
   end

   %- interpolate topo at node locations
   topo_nodes(:,2)=interp1(topo(:,1),topo(:,2),topo_nodes(:,1),'linear');
end

% store in mesh structure
mesh.topo_nodes=topo_nodes;

if(input.debug==1) save('topo_nodes.mat','topo_nodes');end


%-- correct z-coordinates of electrodes for topography
for i=1:length(topo_nodes(:,1))
    for j=1:length(mesh.orig_probe_x)
        if abs(mesh.orig_probe_x(j)-topo_nodes(i,1)) < input.accuracy_threshold
           mesh.orig_probe_z(j)=mesh.orig_probe_z(j)-topo_nodes(i,2);
        end
    end
end


%-- correct z-coordinates of nodes for topography (ie translate ALL nodes)
tmp=unique(node(:,2));
z_except1=tmp(end-1);   % bottom of physical model
z_except2=bottom;       % bottom of total model, including boundary conditions
for i=1:length(node(:,1))
    flag_translated=0;
    for j=1:length(topo_nodes)
        if node(i,1)==topo_nodes(j,1) && node(i,2)==z_except1
        %- translate bottom of physical model by constant shift (min. topo)
           node(i,2)=node(i,2)-min(topo_nodes(:,2));
           flag_translated=1;
        elseif node(i,1)==topo_nodes(j,1) && node(i,2)==z_except2
        %- translate bottom of total model by constant shift (mean topo)
           node(i,2)=node(i,2)-mean(topo_nodes(:,2));
           flag_translated=1;
        elseif node(i,1)==topo_nodes(j,1)
        %- translate other nodes by topo shift
           node(i,2)=node(i,2)-topo_nodes(j,2);
           flag_translated=1;
        end
    end

    %- check if node has well been translated
    if input.topo_flag==1 && flag_translated==0
       disp(sprintf('Node nb %d NOT corrected for topography.', i));
       error('ABORT')
    end
end


%- plot electrode locations
if input.plot_acqui==1
   figure
   plot(mesh.orig_probe_x,-mesh.orig_probe_z,'o');
   title('ACQUISITION')
   hold on
end


%- plot locations of all quadripoles
if input.plot_full_acqui==1
   figure
   hold on
   plot(mesh.orig_probe_x,-mesh.orig_probe_z,'o');
   for i=1:input.num_mes
       plot(input.ax(i),-input.az(i)-i,'ro');
       plot(input.bx(i),-input.bz(i)-i,'bo');
       plot(input.mx(i),-input.mz(i)-i,'ko');
       plot(input.nx(i),-input.nz(i)-i,'go');
   end
   title('FULL ACQUISITION')
end


%Here I will create specific element stize on each edge. I am searching for
%edges that are on theleft, right and bottom, and all the rest

%%%h = waitbar(0,'Mesh Design.Please wait up to minutes....');
disp(' ')
disp('Mesh Design. Please wait up to minutes...')
disp(' ')
%
hdata.hmax = input.hmax;
% hdata.hmax = mesh.probe_spacing/3;
% hdata.hmax = input.hstep;
% hdata.edgeh=ll;
% hdata.fun = 
options.output=input.mesh_option_output;
warning off
addpath ('mesh2d');

%-- Generate mesh
if numel(hdata.hmax)==0
   [p,t,fnum]=meshfaces(node,edge,faces,[],options);
else
   [p,t,fnum]=meshfaces(node,edge,faces,hdata,options);
end

%%%waitbar(1/1);
%%%close(h);
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


%-- Find probe coordinates
for i=1:mesh.num_probes

    % init min_dist for Algo 2
    min_dist=100*mesh.probe_spacing;
    min_j=mesh.num_nodes+1;   %init to impossible value
    flag_found=0;

    for j=1:mesh.num_nodes

        %%Algo 1: find node coinciding with electrode location
        %if( (mesh.x_node_coord(j)==mesh.orig_probe_x(i)) && (mesh.y_node_coord(j)==mesh.orig_probe_z(i)))
        %    mesh.mes_nodes(i)=j;
        %end

        %Algo 1b (FL): find node coinciding with electrode location, up to accuracy threshold
        %if( (mesh.x_node_coord(j)==mesh.orig_probe_x(i)) && (mesh.y_node_coord(j)==mesh.orig_probe_z(i)))
        if( abs(mesh.x_node_coord(j)-mesh.orig_probe_x(i))<input.accuracy_threshold && ...
            abs(mesh.y_node_coord(j)-mesh.orig_probe_z(i))<input.accuracy_threshold )
            mesh.mes_nodes(i)=j;
            flag_found=1;
        end

        %%Algo 2 (FL): find closest node to electrode location
        %dist = sqrt( (mesh.x_node_coord(j)-mesh.orig_probe_x(i))^2 + (mesh.y_node_coord(j)-mesh.orig_probe_z(i))^2 );
        %if dist<min_dist
        %   min_dist=dist;
        %   min_j=j;
        %   mesh.mes_nodes(i)=min_j;
        %end

    end   %end for num_nodes

    if flag_found==1
       disp(sprintf('Coordinates found for electrode nb %d, j = %d', i, mesh.mes_nodes(i)));
       %disp(sprintf('Coordinates found for electrode nb %d, j = %d / %d, dist = %f', i, min_j, mesh.mes_nodes(i), min_dist));
    else
       disp(sprintf('Coordinates NOT found for electrode nb %d', i));
       error('ABORT')
    end
end   %end for num_probes
disp(sprintf('\n'))


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


%-- Find EBC nodes
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


%-- Find number of elements in each parameter
for j=1:mesh.num_param
    [la la2]=find(mesh.param(j,:)~=0);
    mesh.no_elm_per_param(j)=length(la2);
end

hold off

end   %end function mesh=create_mesh_gen(input,mesh)
%----------------------------------------------------------------------%

	

%----------------------------------------------------------------------%
%
%   Function CREATE_MESH_GRID: generate arrays of parameter coordinates.
%
%----------------------------------------------------------------------%

function mesh=create_mesh_grid(input,mesh)

%-- Add topography
%load('anaglyfo_data.mat');
topo_nodes=mesh.topo_nodes;
ytopo_param=zeros(mesh.num_param,1);

%-- Define topographic correction for each parameter cell
for i=1:mesh.num_param
    x_tmp=mesh.param_x(i);
    ind=find(topo_nodes(:,1)==x_tmp);
    if length(ind)>0
       ytopo_param(i)=topo_nodes(ind,2);
    else
       error(sprintf('Did not find ytopo for param. i, param_x(i) = %i, %f',i,mesh.param_x(i)));
    end
end


% model limits
mesh.min_x=min(mesh.param_x);
mesh.max_x=max(mesh.param_x);

mesh.min_y=min(mesh.param_y-ytopo_param);
mesh.max_y=max(mesh.param_y-ytopo_param);

%xx=union(mesh.param_x, mesh.param_x); % take the unique x_line
%yy=union(mesh.param_y-ytopo_param, mesh.param_y-ytopo_param); %take the unique y_line
xx=unique(mesh.param_x);       % take the unique x_line
yy=unique(mesh.param_y-ytopo_param); % take the unique y_line and correct for topo

dx=xx(2)-xx(1);   % find the x step
dy=yy(2)-yy(1);   % find the (min.) y_step
yy=-yy;           % create y as negative

% create grid data
[mesh.XI mesh.YI]=meshgrid(xx,yy);

xnew=[mesh.min_x:dx/2:mesh.max_x];
ynew=[mesh.min_y:dy/2:mesh.max_y];

ynew=-ynew; %Anapoda giati einai arnitika
[mesh.xxx,mesh.yyy]=meshgrid(xnew,ynew);

% store topographic correction for each parameter (useful for later plots)
mesh.param_ytopo=ytopo_param;


%-- Define neighbours
%FL, Sept. 2016: actually not used yet, but potentially useful...
yy=unique(mesh.param_y); % take the unique y_line (without topographic correction)
nx=length(xx);
ny=length(yy);
% loop over parameter cells
for i=1:mesh.num_param

  % unique indices of cell
  ix=find( xx==mesh.param_x(i) );
  iy=find( yy==mesh.param_y(i) );

  % define ranges of x and y where to search for neighbours
  % (taking boundaries into account)
  xmin_search=xx( max(1, ix-1) );
  xmax_search=xx( min(nx,ix+1) );
  ymin_search=yy( max(1, iy-1) );
  ymax_search=yy( min(ny,iy+1) );

  % indices of neighbours
  i_nbrs=find( mesh.param_x>=xmin_search & mesh.param_x<=xmax_search ...
             & mesh.param_y>=ymin_search & mesh.param_y<=ymax_search ...
             & ( mesh.param_x~=xx(ix) | mesh.param_y~=yy(iy) )       );
             % (do not include cell in the list of neighbours)
  mesh.nbrs{i}=i_nbrs;
end


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

end   %-- end function mesh=create_mesh_grid(input,mesh)
%----------------------------------------------------------------------%



%----------------------------------------------------------------------%
%
%   FUNCTION geometrical_factor: compute the geometrical factor K for
%                                each quadripole of electrodes.
%
%----------------------------------------------------------------------%

function mesh=geometrical_factor(input,mesh)

% include first topography in ax and az coordinates
% for i=1:input.num_mes
%     for j=1:length(mesh.topo_nodes(:,1))
%         if input.ax(i)==mesh.topo_nodes(j,1)
%             input.az(i)=input.az(i)-topo_nodes(j,2);
%         end
%         if input.bx(i)==mesh.topo_nodes(j,1)
%             input.bz(i)=input.bz(i)-topo_nodes(j,2);
%         end
%         if input.mx(i)==mesh.topo_nodes(j,1)
%             input.mz(i)=input.mz(i)-topo_nodes(j,2);
%         end
%         if input.nx(i)==mesh.topo_nodes(j,1)
%             input.nz(i)=input.nz(i)-topo_nodes(j,2);
%         end
%     end
% end

% init.
VAM=zeros(input.num_mes,1);
VAN=zeros(input.num_mes,1);
VBM=zeros(input.num_mes,1);
VBN=zeros(input.num_mes,1);


%-- loop over data points/quadripoles
for i=1:input.num_mes
    if(input.m_array_type(i)==1 || input.m_array_type(i)==2)
         aflag=1; bflag=1; mflag=1; nflag=1;
    end

    %- flags for Pole-Dipole array
    if(input.m_array_type(i)==2)
        aflag=1; bflag=0; mflag=1; nflag=1;
    end

    %- flags for Pole-Pole array
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
end   %-- end loop over data points/quadripoles


%-- compute geometrical factor
mesh.geofac=1./(VAM-VBM-VAN+VBN);
clear VAM;clear VBM;clear VAN;clear VBN;clear r1;clear r2;

if input.debug==1
%- output geometrical factor for check
   Mts=[[1:length(mesh.geofac)]' mesh.geofac];
   save('geofac.dat','Mts','-ascii')
end


if input.data_type==2
%if input data are resistances, convert them to resistivities
   input.real_data=input.real_data.*mesh.geofac;
   %pause
end


if input.time_lapse_flag==1
    ind=find (input.d4_real_data<0);
    if size(ind,1)~=0 && size(ind,2)~=-0
        display('Negative resistivities in measurements');
        ind
        pause
    end

%FL /!\ Default values must match homogeneous background,
%       otherwise solution is wrong
% (now defined in separate routine define_mean_res.m)
%%    mesh.mean_res=mean(mean(10.^(log10(input.d4_real_data))));  

elseif input.time_lapse_flag==0
    ind=find (input.real_data<0);
    if size(ind,1)~=0 && size(ind,2)~=-0
        display('Negative resistivities in measurements');
        ind
        pause
    end

%FL /!\ Default values must match homogeneous background,
%       otherwise solution is wrong
% (now defined in separate routine define_mean_res.m)
%%    mesh.mean_res=(mean(10.^(log10(input.real_data))));

end   %end if time_lapse_flag

end   %-- end function mesh=geometrical_factor(input,mesh)
%----------------------------------------------------------------------%

