function [mesh,input]=array3(input)
% global input d4
%==========================================================================
%                  FIND PROBE_X PROBE_Y PROBE_Z                           %
%==========================================================================
% find all coordinates
all_electrodes=[input.ax input.ay input.az];
all_electrodes=[all_electrodes; input.bx input.by input.bz];
all_electrodes=[all_electrodes; input.mx input.my input.mz];
all_electrodes=[all_electrodes; input.nx input.ny input.nz];
all_electrodes=union(all_electrodes,all_electrodes,'rows');

mesh.probe_x=all_electrodes(:,1);
mesh.probe_y=all_electrodes(:,2);
mesh.probe_z=all_electrodes(:,3);

mesh.min_probe_x=min(mesh.probe_x);
mesh.min_probe_y=min(mesh.probe_y);
mesh.min_probe_z=min(mesh.probe_z);
%Save Electrode Coordinates
mesh.X_probe=mesh.probe_x;
mesh.Y_probe=mesh.probe_y;
mesh.Z_probe=mesh.probe_z;

% figure
% plot3(mesh.probe_x,mesh.probe_y,mesh.probe_z,'o');
% Find tolal number of probes
mesh.num_probes=length(mesh.probe_x);

clear all_electrodes



%--------------------------------------------------------------------------
% Find probe_spacing_x and number of different electrodes in x-direction
x_el=unique(mesh.probe_x);
mesh.num_x_probes=length(x_el);
tmp_spacing=zeros(mesh.num_x_probes-1,1);
for i=1:mesh.num_x_probes-1
    tmp_spacing(i)=x_el(i+1)-x_el(i);
end
mesh.probe_spacing_x=min(tmp_spacing);           
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Find probe_spacing_y and number of different electrodes in y-direction
y_el=unique(mesh.probe_y);
mesh.num_y_probes=length(y_el);
tmp_spacing=zeros(mesh.num_y_probes-1,1);
for i=1:mesh.num_y_probes-1
    tmp_spacing(i)=y_el(i+1)-y_el(i);
end
mesh.probe_spacing_y=min(tmp_spacing);       
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Find probe_spacing_z and number of different electrodes in z-direction
z_el=unique(mesh.probe_z);
mesh.num_z_probes=length(z_el);
% Check Here if Only surface electrodes are found
if mesh.num_z_probes==1 %only 0 as coordinate==>Surface
        mesh.num_z_probes=0;
        mesh.probe_spacing_z=0;
        mesh.num_br_probes=0;
else    
        tmp_spacing=zeros(mesh.num_z_probes-1,1);
    for i=1:mesh.num_z_probes-1
        tmp_spacing(i)=z_el(i+1)-z_el(i);
    end
    mesh.probe_spacing_z=min(tmp_spacing);
    
    % Find Number of Borehole Probes (exclude the z=0)
    mesh.num_br_probes=length(z_el);
    if min(z_el)==0 ;mesh.num_br_probes=mesh.num_br_probes-1; end
end
%--------------------------------------------------------------------------
clear tmp_spacing all_electrodes x_el y_el z_el




%==========================================================================
%     assign initial dummy values to probe coord and probe positions      %
%==========================================================================
mesh.pa=zeros(input.num_mes,1);
mesh.pb=zeros(input.num_mes,1);
mesh.pm=zeros(input.num_mes,1);
mesh.pn=zeros(input.num_mes,1);
%==========================================================================
%            Electrode Location for Each Measurement                      %
%==========================================================================
for i=1:input.num_mes
    for j=1:mesh.num_probes	  
        if((input.ax(i)==mesh.probe_x(j)) && (input.ay(i)==mesh.probe_y(j)) && (input.az(i)==mesh.probe_z(j))) ;mesh.pa(i)=int32(j);end
        if((input.bx(i)==mesh.probe_x(j)) && (input.by(i)==mesh.probe_y(j)) && (input.bz(i)==mesh.probe_z(j))) ;mesh.pb(i)=int32(j);end
        if((input.mx(i)==mesh.probe_x(j)) && (input.my(i)==mesh.probe_y(j)) && (input.mz(i)==mesh.probe_z(j))) ;mesh.pm(i)=int32(j);end
        if((input.nx(i)==mesh.probe_x(j)) && (input.ny(i)==mesh.probe_y(j)) && (input.nz(i)==mesh.probe_z(j))) ;mesh.pn(i)=int32(j);end
        if(input.m_array_type(i)==2 && input.bx(i)==0 && input.by(i)==0 && input.bz(i)==0) ;mesh.pb(i)=0;end
%         if(input.m_array_type(i)==3 && input.nx(i)==0 && input.ny(i)==0 && input.nz(i)==0) ;mesh.pn(i)=0;end
        if(input.m_array_type(i)==3) ;mesh.pb(i)=0; mesh.pn(i)=0; end
    end
end


%==========================================================================
%                   GENERALIZED GEOMETRIC FACTOR
%==========================================================================
mesh.geofac=zeros(input.num_mes,1);
 for i=1:input.num_mes
   
	  VAM=0; VAN=0; VBM=0; VBN=0;
      if(mesh.pa(i)*mesh.pm(i)~=0)
			 r1=sqrt( (input.ax(i)-input.mx(i))*(input.ax(i)-input.mx(i))+ ...
				  (input.ay(i)-input.my(i))*(input.ay(i)-input.my(i))+ ...
				  (input.az(i)-input.mz(i))*(input.az(i)-input.mz(i)));

			 r2=sqrt( (input.ax(i)-input.mx(i))*(input.ax(i)-input.mx(i))+ ...
				  (input.ay(i)-input.my(i))*(input.ay(i)-input.my(i))+ ...
				  (-input.az(i)-input.mz(i))*(-input.az(i)-input.mz(i)));

			 VAM=(r1+r2)/(r1*r2*4*pi);
      end

      if(mesh.pa(i)*mesh.pn(i)~=0)
			 r1=sqrt( (input.ax(i)-input.nx(i))*(input.ax(i)-input.nx(i))+ ...
				  (input.ay(i)-input.ny(i))*(input.ay(i)-input.ny(i))+ ...
				  (input.az(i)-input.nz(i))*(input.az(i)-input.nz(i)));

			 r2=sqrt( (input.ax(i)-input.nx(i))*(input.ax(i)-input.nx(i))+ ...
				  (input.ay(i)-input.ny(i))*(input.ay(i)-input.ny(i))+ ...
				  (-input.az(i)-input.nz(i))*(-input.az(i)-input.nz(i)));

			 VAN=(r1+r2)/(r1*r2*4*pi);
      end

	  if (mesh.pb(i)*mesh.pm(i)~=0)
			 r1=sqrt( (input.bx(i)-input.mx(i))*(input.bx(i)-input.mx(i))+ ...
				  (input.by(i)-input.my(i))*(input.by(i)-input.my(i))+ ...
				  (input.bz(i)-input.mz(i))*(input.bz(i)-input.mz(i)));

			 r2=sqrt( (input.bx(i)-input.mx(i))*(input.bx(i)-input.mx(i))+ ...
				  (input.by(i)-input.my(i))*(input.by(i)-input.my(i))+ ...
				  (-input.bz(i)-input.mz(i))*(-input.bz(i)-input.mz(i)));

			 VBM=(r1+r2)/(r1*r2*4*pi);
      end

      if(mesh.pb(i)*mesh.pn(i)~=0)
			 r1=sqrt( (input.bx(i)-input.nx(i))*(input.bx(i)-input.nx(i))+ ...
				  (input.by(i)-input.ny(i))*(input.by(i)-input.ny(i))+ ...
				  (input.bz(i)-input.nz(i))*(input.bz(i)-input.nz(i)));

			 r2=sqrt( (input.bx(i)-input.nx(i))*(input.bx(i)-input.nx(i))+ ...
				   (input.by(i)-input.ny(i))*(input.by(i)-input.ny(i))+ ...
				  (-input.bz(i)-input.nz(i))*(-input.bz(i)-input.nz(i)));
	
			 VBN=(r1+r2)/(r1*r2*4*pi);
      end
      mesh.geofac(i)=1/(VAM-VBM-VAN+VBN);
 
 end
%==================================================================================================
% find ap. res if initial data given in resistances */
      
if (input.data_type==2)
      input.real_data=input.real_data.*mesh.geofac;       
end
if input.time_lapse_flag==0; input.old_real_data=input.real_data; end

 %=========================================================================
 %  Fine initial model
 if input.time_lapse_flag==0
 mesh.mean_res=10^mean(log10(input.real_data));
 elseif input.time_lapse_flag==1
     mesh.mean_res=10.^mean(log10(input.d4_real_data));
     mesh.mean_res=mean(mesh.mean_res);
 end
%    mesh.mean_res=20;
% mesh.mean_res=complex(10,eps);
% input.sip_flag=1;
%==================================================================================================


%===========================================================================
%            Maximum Nsep and Number of Parameter Layers
%===========================================================================
% find max_n or maximum layers */
if (mesh.num_z_probes~=0)
    mesh.n_sep=mesh.num_z_probes; 
elseif mesh.num_z_probes==0
max_sr_n=max(sqrt((input.ax-input.mx).^2+(input.az-input.mz).^2));
max_sr_n=int32(max_sr_n/   (min (mesh.probe_spacing_x,mesh.probe_spacing_y)));     
  %N-SEP value for the optimized data sets
%   mesh.n_sep=mesh.num_x_probes-1;
  mesh.n_sep=6;
  mesh.max_n=mesh.n_sep;

end



end