function mesh=mesh_gen_3d_new_version(input,mesh,depth_n)
% global input cm d4


num_of_hex=2; %or 4

num_ext=14;
num_bot=14;
right=[1 1 1 2 3 5 10 25 32 50 100 250 400 800];
% left= [0 25 10 5 3 2 1];
bottom=[1 1 1 2 3 5 10 25 32 50 100 250 400 800];
right=right.*mesh.probe_spacing_x/2;%*100;
% left=left.*mesh.probe_spacing_x;%*100;
bottom=bottom.*mesh.probe_spacing_x/2;%*100;

%/* find the min,mid,max dimensions */
mesh.scale_x=mesh.probe_spacing_x/mesh.probe_spacing_y;
mesh.scale_y=mesh.probe_spacing_y/mesh.probe_spacing_x;
mesh.scale_z=mesh.probe_spacing_z/mesh.probe_spacing_x;
if(mesh.scale_z==0)  ;mesh.scale_z=0.5; end 




% The main idea is to find the unique x,y,x points based on those I could
% create the mesh. I need first to find the UNIQUE electrodes coordinates (see
% function array3) so i have unique_x unique_y unique_z.
% Now in each of those vector I need to extent them (boundary conitions)
% and create one mid point between two electrode positions.




%* *******  generate base coordinates ************** */
     %--------------------------------------------------------------------%
     %/* base x */
     clear base_x tmp_matrix
     unique_x=unique(mesh.probe_x);     
    
     % interpote mid points
     %unique_x=[0:0.25:12];
     
     
     
     reference_start=min(unique_x); % From this point start REDUCING
     base_x(num_ext+1)=reference_start;
     for i=1:num_ext
         base_x(num_ext+1-i)=base_x(num_ext+2-i)-right(i);
     end
     
     left_x=unique_x(1);
     tmp1=1;
     for i=1:length(unique_x)-1
         tmp_matrix(tmp1)=unique_x(i);
         tmp1=tmp1+1;
         tmp_matrix(tmp1)=(unique_x(i)+unique_x(i+1) )/2;
         tmp1=tmp1+1;
     end
     tmp_matrix(tmp1)=max(unique_x);
     
     base_x=union(tmp_matrix,base_x);
  	
     tmp=1;
     for i=length(base_x)+1:length(base_x)+num_ext
        base_x(i)=base_x(i-1)+right(tmp);
        tmp=tmp+1;
     end
     tmpx=length(base_x);
	 mesh.xc_mesh=(base_x(1)+base_x(tmpx))/2.;
     %--------------------------------------------------------------------%
     %/* base y */
     unique_y=unique(mesh.probe_y);
     %unique_y=[0:0.25:12];
     
     clear base_y tmp_matrix
    
     
     reference_start=min(unique_y); % From this point start REDUCING
     base_y(num_ext+1)=reference_start;
     for i=1:num_ext
         base_y(num_ext+1-i)=base_y(num_ext+2-i)-right(i)*mesh.scale_y;
     end
     
     left_y=unique_y(1);
     tmp1=1;
     for i=1:length(unique_y)-1
         tmp_matrix(tmp1)=unique_y(i);
         tmp1=tmp1+1;
         tmp_matrix(tmp1)=(unique_y(i)+unique_y(i+1) )/2;
         tmp1=tmp1+1;
     end
     tmp_matrix(tmp1)=max(unique_y);
     
     base_y=union(tmp_matrix,base_y);
  	
     tmp=1;
     for i=length(base_y)+1:length(base_y)+num_ext
        base_y(i)=base_y(i-1)+right(tmp);
        tmp=tmp+1;
     end
     tmpy=length(base_y);
	 mesh.yc_mesh=(base_y(1)+base_y(tmpy))/2.;
     
%--------------------------------------------------------------------------     
     %/*base z*/
	 clear base_z
     base_z(1)=0;
	 
	 tmp=0;
	unique_z=unique(mesh.probe_z);
    unique_z(1)=0;
      
    
     if mesh.num_z_probes==0
        probe_spacing=min(mesh.probe_spacing_x,mesh.probe_spacing_y);
         unique_z(2)=probe_spacing/2;
         m_factor=1.15; % This means 15%thicker each layer
         for j=3:mesh.n_sep+1 
             tmp_thick=unique_z(j-1)-unique_z(j-2);
             unique_z(j)=unique_z(j-1)+m_factor*tmp_thick;
         end
         
         % add here info for custom depth_n
         if length(depth_n)>1
            unique_z=depth_n';
            mesh.n_sep=length(depth_n)-1;
         end
         
         
         %unique_z=[0:0.125:4];
         %depth_n=unique_z;
         %mesh.n_sep=length(depth_n)-1;
      
      
      
         tmp1=1;
         for j=1:mesh.n_sep
             base_z(tmp1)=unique_z(j);
             tmp1=tmp1+1;
             base_z(tmp1)=(unique_z(j)+unique_z(j+1))/2;
             tmp1=tmp1+1;
         end
         base_z(tmp1)=unique_z(length(unique_z));
%          for j=2:mesh.n_sep+1 ;unique_z(j)=unique_z(j-1)+mesh.probe_spacing_y/2;end

%         for i=2:(2*mesh.n_sep+1) ;base_z(i)=base_z(i-1)+mesh.probe_spacing_y/4;end
        
        
        
     else
         unique_z=union(unique_z,0); % make sure we have zero coordinate
%          for j=1:mesh.n_sep+1 ;base_probe_z(j)=unique_z(j);end
         tmp1=1;
         for i=1:length(unique_z)-1
             base_z(tmp1)=(unique_z(i));
             tmp1=tmp1+1;
             base_z(tmp1)=( unique_z(i)+unique_z(i+1) )/2;
             tmp1=tmp1+1;
         end
         base_z(tmp1)=max(unique_z);
     end
     
     for i=2*mesh.n_sep+2:2*mesh.n_sep+num_bot+1
        tmp=tmp+1;
        base_z(i)=base_z(i-1)+bottom(tmp)*mesh.scale_z;
     end
     tmpz=length(base_z);

maxbz=base_z(i);

mesh.depth_n=unique_z';
tmp=unique_z';
save('layers.mat','tmp');

%--------------------------------------------------------------------------
%                                MESHGRID STUF
%--------------------------------------------------------------------------
% I could replace a lot of stuff by using
%[XI YI ZI]=meshgrod(base_x,base_y,base_z).
% Maybe in future.


tmp=max(tmpx,tmpy); max_d=max(tmp,tmpz);
tmp=min(tmpx,tmpy); min_d=min(tmp,tmpz);

     if((tmpx==min_d) && (tmpy==max_d)) ;mid_d=tmpz; fz=2;end
     if((tmpx==min_d) && (tmpz==max_d)) ;mid_d=tmpy; fy=2;end
     if((tmpy==min_d) && (tmpz==max_d)) ;mid_d=tmpx; fx=2;end

     if((tmpx==max_d) && (tmpy==min_d)) ;mid_d=tmpz; fz=2;end
     if((tmpx==max_d) && (tmpz==min_d)) ;mid_d=tmpy; fy=2;end
     if((tmpy==max_d) && (tmpz==min_d)) ;mid_d=tmpx; fx=2;end

     if(max_d==tmpx) ;fx=3;end; if(min_d==tmpx) ;fx=1;end
     if(max_d==tmpy) ;fy=3;end; if(min_d==tmpy) ;fy=1;end
     if(max_d==tmpz) ;fz=3;end; if(min_d==tmpz) ;fz=1;end
     
       %/* check for possible equal flags */
      if((fx==fy) && (fx==fz))  ;fx=3; max_d=tmpx;  fy=2; mid_d=tmpy; fz=1; min_d=tmpz;end
      if((fx==fy) && (fx>fz))   ;fx=3; max_d=tmpx;  fy=2; mid_d=tmpy;end
      if((fx==fy) && (fx<fz))   ;fx=2; mid_d=tmpx;  fy=1; min_d=tmpy;end
      if((fz==fy) && (fz>fx))   ;fy=3; max_d=tmpy;  fz=2; mid_d=tmpz;end
      if((fz==fy) && (fz<fx))   ;fy=2; mid_d=tmpy;  fz=1; min_d=tmpz;end
      if((fx==fz) && (fx>fy))   ;fx=3; max_d=tmpx;  fz=2; mid_d=tmpz;end
      if((fx==fz) && (fx<fy))   ;fx=2; mid_d=tmpx;  fz=1; min_d=tmpz;end
      
      
% Really important. Ensurance pc accurancy
base_y=round(base_y*1000)/1000;
base_x=round(base_x*1000)/1000;
base_z=round(base_z*1000)/1000;
mesh.x_node_coord=zeros(max_d*mid_d*min_d,1);
mesh.y_node_coord=zeros(max_d*mid_d*min_d,1);
mesh.z_node_coord=zeros(max_d*mid_d*min_d,1);
       %/* ***** find nodal coordinates ****** */
	   l=0;
	   for i=1:max_d
	      for k=1:mid_d
		     for j=1:min_d
		  
		l=l+1;
		if(fx==1) ;mesh.x_node_coord(l)=base_x(j); end
		if(fy==1) ;mesh.y_node_coord(l)=base_y(j); end
		if(fz==1) ;mesh.z_node_coord(l)=base_z(j); end

		if(fx==2) ;mesh.x_node_coord(l)=base_x(k); end
		if(fy==2) ;mesh.y_node_coord(l)=base_y(k); end
		if(fz==2) ;mesh.z_node_coord(l)=base_z(k); end

		if(fx==3) ;mesh.x_node_coord(l)=base_x(i); end
		if(fy==3) ;mesh.y_node_coord(l)=base_y(i); end
		if(fz==3) ;mesh.z_node_coord(l)=base_z(i); end
             end
          end
       end
	       mesh.num_nodes=l;
%--------------------------------------------------------------------------
%              %* ********** find measurement nodes*************** */
%--------------------------------------------------------------------------
mesh.mes_nodes=[];
%   l=0;
%   for j=1:mesh.num_probes
%    for i=1:mesh.num_nodes
%      if( (mesh.probe_x(j)==mesh.x_node_coord(i)) && (mesh.probe_y(j)==mesh.y_node_coord(i)) && (mesh.probe_z(j)==mesh.z_node_coord(i)))
%       l=l+1;
%       mesh.mes_nodes(l)=i;
% %       mesh.mes_nodes(l,2)=j;
%      end
%    end
%   end
  
  
  for j=1:mesh.num_probes
     ind=find(mesh.probe_x(j)==mesh.x_node_coord & mesh.probe_y(j)==mesh.y_node_coord & mesh.probe_z(j)==mesh.z_node_coord);
     mesh.mes_nodes(j)=ind;
      
  end
  
  
  
if length(mesh.mes_nodes)~=mesh.num_probes
    disp('Something went wrong');
    pause
end

  %* ********** find ebc nodes *************** */

% locate initally boarder nodes and store them tempor. in icon[][] */
% mesh.node_ebc=[];
%     k=0;
% 
%     for i=1:mesh.num_nodes
%        
% 	 l=0;
% 	 if(mesh.x_node_coord(i)==base_x(1))     ;k=k+1; l=l+1; mesh.node_ebc(k)=i; end
% 	 if(mesh.x_node_coord(i)==base_x(tmpx))  ;k=k+1; l=l+1; mesh.node_ebc(k)=i; end
% 
% 	 if(mesh.y_node_coord(i)==base_y(1))	    ;k=k+1; l=l+1; mesh.node_ebc(k)=i; end
% 	 if(mesh.y_node_coord(i)==base_y(tmpy))  ;k=k+1; l=l+1; mesh.node_ebc(k)=i; end
% 
% 	 if(mesh.z_node_coord(i)==base_z(tmpz))  ;k=k+1; l=l+1; mesh.node_ebc(k)=i; end
% 
% 	 if(l~=0) ;k=k-l+1; end
%     end
%       mesh.num_ebc=k;
%       mesh.node_ebc=mesh.node_ebc(1:k);
      
      
      ind1=find(mesh.x_node_coord==base_x(1));
      ind2=find(mesh.x_node_coord==base_x(tmpx));
      ind3=find(mesh.y_node_coord==base_y(1));
      ind4=find(mesh.y_node_coord==base_y(tmpy));
      ind5=find(mesh.z_node_coord==base_z(tmpz));
      
      ind=union(ind1,ind2);
      ind=union(ind,ind3);
      ind=union(ind,ind4);
      ind=union(ind,ind5);
      ind=unique(ind);
      
      mesh.node_ebc=ind;
      mesh.num_ebc=length(ind);
clear ind ind1 ind2 ind3 ind4 ind5      
      
  %* ******** generate mesh information *************** */
mesh.icon=[];
   cnt_icon=0;
   for l=1:max_d-1
    k=(l-1)*min_d*mid_d;
	for i=1:mid_d-1	   
	   for j=k+(i-1)*min_d+1:k+i*min_d-1
	   		cnt_icon=cnt_icon+1;
            mesh.icon(1,cnt_icon)=j;
            mesh.icon(2,cnt_icon)=j+min_d;
            mesh.icon(3,cnt_icon)=j+min_d*mid_d;
            mesh.icon(4,cnt_icon)=j+(min_d*mid_d)+min_d;
            mesh.icon(5,cnt_icon)=j+1;
            mesh.icon(6,cnt_icon)=j+min_d+1;
            mesh.icon(7,cnt_icon)=j+min_d*mid_d+1;
            mesh.icon(8,cnt_icon)=j+min_d*mid_d+min_d+1;
            mesh.prop(cnt_icon)=mesh.mean_res;
        end
    end
   end
	      mesh.num_elements=cnt_icon;

% /* ****************** define parameters ************************ */
% /*  param_i(elem)(1)=x coord of param, (2)=y coor (3)=z coor
%      (4)(5) up down x_lim, (6)(7) up down y_lim (8)(9) up down z_lim*/
% /* ************************************************************ */  

%  /* define parameter mesh coordinates and search limits */  
mesh.param_i=[];
mesh.param_x=[];
mesh.param_y=[];
mesh.param_z=[];
disp('     DEFINE 3D PARAMETERS MATRIX..');tic
  m=0;
  for i=1:mesh.n_sep
    for j=1:length(unique_y)-1
	  for k=1:length(unique_x)-1

		 m=m+1;
		 mesh.param_i(m,1)=double((unique_x(k)+unique_x(k+1))/2);
		 mesh.param_i(m,2)=(unique_y(j)+unique_y(j+1))/2;
		 mesh.param_i(m,3)=(unique_z(i)+unique_z(i+1))/2;

		 mesh.param_x(m,1)=(mesh.param_i(m,1)-left_x);
		 mesh.param_y(m,1)=(mesh.param_i(m,2)-left_y);
		 mesh.param_z(m,1)=(mesh.param_i(m,3));


		 mesh.param_i(m,4)=unique_x(k);
		 mesh.param_i(m,5)=unique_x(k+1);
		 mesh.param_i(m,6)=unique_y(j);
		 mesh.param_i(m,7)=unique_y(j+1);
		 mesh.param_i(m,8)=unique_z(i);
		 mesh.param_i(m,9)=unique_z(i+1);

		%if(opt_flag==0)
		
		  if(k==1) ;mesh.param_i(m,4)=min(base_x); end
		  if(k==length(unique_x)-1) ;mesh.param_i(m,5)=max(base_x); end
		  if(j==1) ;mesh.param_i(m,6)=min(base_y); end
		  if(j==(length(unique_y)-1)) ;mesh.param_i(m,7)=max(base_y); end
		  if(i==1) ;mesh.param_i(m,8)=0; end
		  if(i==mesh.n_sep) ;mesh.param_i(m,9)=max(base_z); end
        %end
      end
    end
  end
	     mesh.num_param=m;

% /* assign initial resistivity to parameters */
mesh.res_param1(1:mesh.num_param,1)=mesh.mean_res;
mesh.res_param2(1:mesh.num_param,1)=mesh.mean_res;
if input.time_lapse_flag==1
   mesh.d4_res_param1(1:mesh.num_param,1:input.num_files)=mesh.mean_res;
   mesh.d4_res_param2=mesh.d4_res_param1;    
end





mesh.elem_param=zeros(mesh.num_elements,1);
mesh.jtmp=zeros(1500,mesh.num_param);
% num_of_elements=zeros(mesh.num_param,1);


	   e1=mesh.icon(1,:); e2=mesh.icon(2,:); e3=mesh.icon(3,:); e4=mesh.icon(4,:);
	   e5=mesh.icon(5,:); e6=mesh.icon(6,:); e7=mesh.icon(7,:); e8=mesh.icon(8,:);

       
  xe=(mesh.x_node_coord(e1)+mesh.x_node_coord(e2)+mesh.x_node_coord(e3)+mesh.x_node_coord(e4)+ ...
      mesh.x_node_coord(e5)+mesh.x_node_coord(e6)+mesh.x_node_coord(e7)+mesh.x_node_coord(e8))/8;
  ye=(mesh.y_node_coord(e1)+mesh.y_node_coord(e2)+mesh.y_node_coord(e3)+mesh.y_node_coord(e4)+ ...
      mesh.y_node_coord(e5)+mesh.y_node_coord(e6)+mesh.y_node_coord(e7)+mesh.y_node_coord(e8))/8;
  ze=(mesh.z_node_coord(e1)+mesh.z_node_coord(e2)+mesh.z_node_coord(e3)+mesh.z_node_coord(e4)+ ...
      mesh.z_node_coord(e5)+mesh.z_node_coord(e6)+mesh.z_node_coord(e7)+mesh.z_node_coord(e8))/8;
       
    for j=1:mesh.num_param
                            
             ind=find((xe>=mesh.param_i(j,4)) & (xe<=mesh.param_i(j,5)) & ...
                 (ye>=mesh.param_i(j,6)) & (ye<=mesh.param_i(j,7)) & ...
                     (ze>=mesh.param_i(j,8)) & (ze<=mesh.param_i(j,9))) ;
				   
             mesh.elem_param(ind)=j; %conduct search to locate the elements of each parameter */
             mesh.param_i(j,10)=length(ind); % Find here which elements belong to each parameter
             mesh.jtmp(1:length(ind),j)=ind;
             
    end

       

    

    
    
%conduct search to locate the elements of each parameter */
% mesh.elem_param2=zeros(mesh.num_elements,1);
%       for i=1:mesh.num_elements
% 	   % find x y z centr of elem */
% 	   e1=mesh.icon(1,i); e2=mesh.icon(2,i); e3=mesh.icon(3,i); e4=mesh.icon(4,i);
% 	   e5=mesh.icon(5,i); e6=mesh.icon(6,i); e7=mesh.icon(7,i); e8=mesh.icon(8,i);
% 
%   xe=(mesh.x_node_coord(e1)+mesh.x_node_coord(e2)+mesh.x_node_coord(e3)+mesh.x_node_coord(e4)+ ...
%       mesh.x_node_coord(e5)+mesh.x_node_coord(e6)+mesh.x_node_coord(e7)+mesh.x_node_coord(e8))/8;
%   ye=(mesh.y_node_coord(e1)+mesh.y_node_coord(e2)+mesh.y_node_coord(e3)+mesh.y_node_coord(e4)+ ...
%       mesh.y_node_coord(e5)+mesh.y_node_coord(e6)+mesh.y_node_coord(e7)+mesh.y_node_coord(e8))/8;
%   ze=(mesh.z_node_coord(e1)+mesh.z_node_coord(e2)+mesh.z_node_coord(e3)+mesh.z_node_coord(e4)+ ...
%       mesh.z_node_coord(e5)+mesh.z_node_coord(e6)+mesh.z_node_coord(e7)+mesh.z_node_coord(e8))/8;
% 
% 
% 	 for j=1:mesh.num_param
% 	   xf=0; yf=0; zf=0; 
% % % 	 /*  if((mesh.param_i(j,4)<=xe) && (mesh.param_i(j,5)>=xe)) xf=1;
% % % 	   if((mesh.param_i(j,6)<=ye) && (mesh.param_i(j,7)>=ye)) yf=1;
% % % 	   if((mesh.param_i(j,8)<=ze) && (mesh.param_i(j,9)>=ze)) zf=1;*/
% 
% 	if((xe>=mesh.param_i(j,4)) && (xe<=mesh.param_i(j,5)) && ...
% 	   (ye>=mesh.param_i(j,6)) && (ye<=mesh.param_i(j,7)) && ...
% 	   (ze>=mesh.param_i(j,8)) && (ze<=mesh.param_i(j,9))) 
% 				    mesh.elem_param2(i)=j;  break;
%     end
% %    if xf*yf*zf==1 ;mesh.elem_param(i)=j; break; end
%      end
%       end
%       
% l=0;
%       for i=1:mesh.num_param
%     	  k=0;
%           for j=1:mesh.num_elements
%               if i==mesh.elem_param(j)
%                   k=k+1; 
%              end
%       end
% 	      mesh.param_i(i,11)=k;
%       end
%   
      
      
      

% Find sink node in center
median_x=round(mean(mesh.x_node_coord)*1000)/1000;
median_y=round(mean(mesh.y_node_coord)*1000)/1000;
mesh.sink=find(mesh.x_node_coord==median_x & mesh.y_node_coord==median_y &mesh.z_node_coord==0);




% % Find here which elements belong to each parameter
% mesh.jtmp2=zeros(1500,mesh.num_param);
% num_of_elements=zeros(mesh.num_param,1);
%     for i=1:mesh.num_param
%         m=0;
%         for k=1:mesh.num_elements
%             if(mesh.elem_param(k)==i) 
%                 m=m+1;
%                 num_of_elements(i)=m;
%                 mesh.jtmp2(m,i)=k;
% 
%             end
%         end
%     end
% % 













%Create model file to use in forward modelling
modout=fopen('modout.mod','wt');
for i=1:mesh.num_param
    fprintf(modout,'%f\t %f\t %f\t 10 \n',mesh.param_x(i),mesh.param_y(i),-mesh.param_z(i));
end
fclose(modout);
end