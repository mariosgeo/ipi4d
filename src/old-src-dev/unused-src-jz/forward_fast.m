function fem=forward_fast(input,mesh,fem,current_file)
% 
% global h2

% matrix preallocation

        fem.array_model_data=zeros(input.num_mes,1);
           fem.normalizing_data=zeros(input.num_mes,1);
        fem.kernel=zeros(3,3,mesh.num_elements);

%         fem.Kp=zeros(mesh.num_nodes,mesh.num_nodes,length(mesh.k));
end

    aaa=zeros(mesh.num_probes,mesh.num_nodes);

tt=sprintf('**  ITERATION =>  %d  **\n',itr);
tt=cellstr(tt);

% if current_file<2
% h2=waitbar(0,'Wave number');
% end


for kit=1:length(mesh.k)

%     if input.time_lapse_flag==0;waitbar(kit/length(mesh.k));end
%     if input.time_lapse_flag==1;waitbar((kit+length(mesh.k)*(current_file-1))/(length(mesh.k)*input.num_files));end
%     drawnow;
    [tmp_aaa,fem]=form(kit,input,mesh,fem,itr);
    fem.tmp_aaa{kit}=tmp_aaa;
    aaa=aaa+mesh.g(kit).*tmp_aaa/pi;% Inverse fourier transform


end
% 
% if input.time_lapse_flag==0; close(h2);end
% if input.time_lapse_flag==1 && ((kit+length(mesh.k)*(current_file-1))==(length(mesh.k)*input.num_files)) ;close(h2);end
% 

fem=save_data(itr,input,mesh,fem,aaa);
end



% 	/* *********************************************************************
% 	   *                              FORM( )                              *
% 	   *-------------------------------------------------------------------*
% 	   * This function forms the global stiffness matrix.                  *
% 	   *                                                                   *
% 	   * This function calculate the stiffnes term for every element.      *
% 	   * It is valid only for linear triangular elements                   *
% 	   * Arguments:                                                        *
% 	   * e= the number of the element                                      *
% 	   *********************************************************************/


function [tmp_aaa,fem]=form(kit,input,mesh,fem,itr)

    a_index=zeros(3*3*mesh.num_elements,1);
    b_index=zeros(3*3*mesh.num_elements,1);
    val_index=zeros(3*3*mesh.num_elements,1);
%     k=(zeros(mesh.num_nodes,mesh.num_nodes));
    f=zeros(mesh.num_nodes,mesh.num_probes);
    

ke=(zeros(3,3));
fe=zeros(3,1);
k_index=1;
% Start to assemble all area elements in Omega 
for e=1:mesh.num_elements

    % Calculate b_node and c_node (i=1,2,3)    
    I=mesh.icon(1,e);
    J=mesh.icon(2,e);
    L=mesh.icon(3,e);
    
   	b_node(1)=mesh.y_node_coord(J)-mesh.y_node_coord(L);
	b_node(2)=mesh.y_node_coord(L)-mesh.y_node_coord(I);
	b_node(3)=mesh.y_node_coord(I)-mesh.y_node_coord(J);

	c_node(1)=mesh.x_node_coord(L)-mesh.x_node_coord(J);
	c_node(2)=mesh.x_node_coord(I)-mesh.x_node_coord(L);
	c_node(3)=mesh.x_node_coord(J)-mesh.x_node_coord(I);
    % Calculate area of element Delta^e
	area_element=(b_node(2)*c_node(3)-b_node(3)*c_node(2))/2;
    % Generate the elemental matrix [K^e]
    for i=1:3
        
%         fe(i)=FINT(elno)*area/3; % look eq 4.34
		for j=i:3
        		if(j==i) factor=6; end
				if(j~=i) factor=12; end
                ke(i,j)=((b_node(i)*b_node(j)+c_node(i)*c_node(j))/(4.*area_element))+...
                                          ((area_element*mesh.k(kit)^2)/(factor));  %look eq 4.33
                                      
%                 me(i,j)=GAMMA(e)*area_element/factor;
                                      
        end
    end
    % Remember that previous double for loop begins with j=i and not j=1
    for i=1:3
        for j=1:3
            ke(j,i)=ke(i,j);
%             me(j,i)=me(i,j);
        end
    end
    
  
   
    
    % Add [K^e] to [K]
    for i=1:3
        ii=mesh.icon(i,e);
%         f(ii)=f(ii)+fe(i);
        for j=1:3
            jj=mesh.icon(j,e);
%             k(mesh.icon(i,e),mesh.icon(j,e))=k(mesh.icon(i,e),mesh.icon(j,e))+ke(i,j)/mesh.prop(e);
            a_index(k_index)=ii;
            b_index(k_index)=jj;
            val_index(k_index)=ke(i,j)/mesh.prop(e);            
            k_index=k_index+1;

        end
    end
    
    
end

k=sparse(a_index,b_index,val_index,mesh.num_nodes,mesh.num_nodes);


% Start with boundary conditions. We have two options, homogeneous
% dirichlet and mixed bounary conditions
mesh.bc_option=1;
if mesh.bc_option==1
    [tmp_aaa,fem]=dirichlet_bc(k,mesh,fem,itr,kit);
elseif mesh.bc_option==2
    [tmp_aaa,fem]=mixed_boundary_conditions(k,mesh,fem,kit,itr);
end




end

% 	/* ********************************************************************
% 	   *                          mixed_boundary_conditions( )            *
% 	   *------------------------------------------------------------------*

function [tmp_aaa,fem]=mixed_boundary_conditions(k,mesh,fem,kit,itr)



k_keep=k; % For each source, I need the original k, before the bc
current=1;
tmp_aaa=zeros(mesh.num_nodes,mesh.num_probes);

for t=1:mesh.num_probes

    f=zeros(mesh.num_nodes,1);
    k=k_keep;

    % find source node

    % current probe location

    % Here assign value for elements that do no have homogeneous Neuman
    % condition (gamma not equal to zero)
    % ns is the number of nodes that do not have homogeneous Neuman
    % ns(1-2,s) something like node_ebc
    % Start to assemble all line segments on Gamma_2 (mixed boundary conditions)
    for s=1:mesh.num_ebc2
           % calculate the length of each segment
           I=mesh.ns(1,s);
           J=mesh.ns(2,s);
           ls=sqrt(  (mesh.x_node_coord(I) - mesh.x_node_coord(J))^2 + (mesh.y_node_coord(I) - mesh.y_node_coord(J))^2);
           % compute [K^s]


           %find x and y middle point
           mid_x=0.5*(mesh.x_node_coord(I)+mesh.x_node_coord(J));
           mid_y=0.5*(mesh.y_node_coord(I)+mesh.y_node_coord(J));
           tmp11=sqrt((mid_x-mesh.orig_probe_x(t)).*(mid_x-mesh.orig_probe_x(t))+ ...
                       (mid_y-mesh.orig_probe_z(t)).*(mid_y-mesh.orig_probe_z(t)) );
           tmp12=sqrt((mid_x-mesh.orig_probe_x(t)).*(mid_x-mesh.orig_probe_x(t))+ ...
                       (-mid_y-mesh.orig_probe_z(t)).*(-mid_y-mesh.orig_probe_z(t)) );



           gamma=(1/mesh.prop(I))*mesh.k(kit)*...
               (  besselk(1,mesh.k(kit).*tmp11) + besselk(1,mesh.k(kit).*tmp12)   ) ./ ...
                             (besselk(0,mesh.k(kit).*tmp11)+besselk(0,mesh.k(kit).*tmp12)  );



           ks(1,1)=gamma*ls/3;
           ks(1,2)=gamma*ls/6;
           ks(2,1)=ks(1,2);
           ks(2,2)=ks(1,1);

           % Add [K^s] to [K]
           for i=1:2
               for j=1:2
                       k(mesh.ns(i,s),mesh.ns(i,j))=k(mesh.ns(i,s),mesh.ns(j,s))+ ks(i,j);
               end
           end



     % Apply the nonzero natural boundary conditions 
     % see equation 4.52 
     % In the resistivity case q=zero, so I don't calculate the loop below.

     % for t=1:mesh.num_probes
         % for n=1:mesh.nub_ebc2
           %     I=ns(1,n);
           %     J=ns(2,n);
           %     ls=sqrt(  (mesh.x_node_coord(I) - mesh.x_node_coord(J))^2 + (mesh.y_node_coord(I) - mesh.y_node_coord(J))^2);
           %     f(I,t)=f(I,t)+q(n)*ls/2;
           %     f(J,t)=f(J,t)+q(n)*ls/2;
         % end
     % end

    end

    
    
    
    
    
    
    
    % keep homogeneous stiffeness matrix, to solve for secondary field.
    if itr<2
%         fem.Kp{kit}=k;
    end       

    % Rememeber that I need to spilt current node to elements.
    % In the DC case (f=a*I/360), we have thow optiops, either node is on
    % surface (a=180), or insede earth where a=360. 
    source=mesh.mes_nodes(t);
    if mesh.y_node_coord(source)~=0
        split=1;
    else
        split=1;
    end
    
    f(source,1)=f(source,1)+current/split;


    
    % regular solution
        tmp_aaa(:,t)=k\f;

    % singularity removal
%     tmp_aaa(:,t)=fem.Ap_bessel(t,:,kit).'+k\((fem.Kp{kit}-k)*fem.Ap_bessel(t,:,kit).');

    % Preconditioned conjugate gradients method

    % opts.michol = 'on';
    % L = ichol(k);
    % [x2,fl2,rr2,it2,rv2] = pcg(k,(fem.Kp(:,:,kit)-k)*fem.Ap_bessel(t,:,kit).',1e-8,100,L,L');
    % [x2,fl2,rr2,it2,rv2] = pcg(k,(fem.Kp(:,:,kit)-k)*fem.Ap_bessel(t,:,kit).',1e-8);

    % tmp_aaa(:,t)=fem.Ap_bessel(t,:,kit).'+cgcomplex(k,(fem.Kp(:,:,kit)-k)*fem.Ap_bessel(t,:,kit).');
    
    
    
    
    
end

tmp_aaa=full(tmp_aaa).';

end

% 	/* ********************************************************************
% 	   *                          homogeneous dirichlet BC( )             *
% 	   *------------------------------------------------------------------*


function [tmp_aaa,fem]=dirichlet_bc(k,mesh,fem,itr,kit)

f=zeros(mesh.num_nodes,mesh.num_probes);

            % Apply the essential boundary conditions (Dirichlet)
% ------------------------------------------------------------------------%

     % This this the matrix that has the non homegeneous dirichlet conditions
     fem.prop_ebc=zeros(mesh.num_ebc,mesh.num_probes);
     % in case that potential is known on BC
%      for t=1:mesh.num_probes
%      tmp11=sqrt((mesh.x_node_coord(mesh.node_ebc)-mesh.orig_probe_x(t)).*(mesh.x_node_coord(mesh.node_ebc)-mesh.orig_probe_x(t))+ ...
%                  (mesh.y_node_coord(mesh.node_ebc)-mesh.orig_probe_z(t)).*(mesh.y_node_coord(mesh.node_ebc)-mesh.orig_probe_z(t)));
%                        
%      tmp12=sqrt((mesh.x_node_coord(mesh.node_ebc)-mesh.orig_probe_x(t)).*(mesh.x_node_coord(mesh.node_ebc)-mesh.orig_probe_x(t))+ ...
%                  (mesh.y_node_coord(mesh.node_ebc)-mesh.orig_probe_z(t)).*(mesh.y_node_coord(mesh.node_ebc)-mesh.orig_probe_z(t))  );
%                        
%      dis=(tmp11.*tmp12./(tmp11+tmp12));
%      fem.prop_ebc(:,t)=mesh.mean_res./(4.*pi) * (besselk(0,mesh.k(kit).*tmp11)+besselk(0,mesh.k(kit).*tmp12)  );
%      end




num_nodes=mesh.num_nodes;
prop_ebc=fem.prop_ebc;
node_ebc=mesh.node_ebc;

for t=1:mesh.num_probes  % do nout use parfor here. running out of memory
                 f(:,t)=f(:,t)-sum( bsxfun(@times,k(:,node_ebc),prop_ebc(:,t).') ,2);
end
       
f(mesh.node_ebc,:)=fem.prop_ebc;
k(mesh.node_ebc,:)=0;
k(:,mesh.node_ebc)=0;
k(sub2ind(size(k),mesh.node_ebc,mesh.node_ebc))=1; 
current=1;


for t=1:mesh.num_probes
    source=mesh.mes_nodes(t);

    % Rememeber that I need to spilt current node to elements.
    % In the DC case (f=a*I/360), we have thow optiops, either node is on
    % surface (a=180), or insede earth where a=360. 
    if mesh.y_node_coord(source)~=0
        split=1;
    else
        split=1;
    end

    f(source,t)=f(source,t)+current/split;
%     f(cm.sink,t)=f(cm.sinc,t)-current;
 end
            
if itr<2
        fem.Kp{kit}=k;
        
        for kit2=1:length(mesh.k)
            for t=1:mesh.num_probes
                fem.Ap_bessel(t,mesh.node_ebc,kit2)=fem.prop_ebc(:,t);
            end
        end
end             
            
%      singularity removal

% for t=1:mesh.num_probes
%     tmp_aaa(:,t)=fem.Ap_bessel(t,:,kit).'+k\((fem.Kp{kit}-k)*fem.Ap_bessel(t,:,kit).');           
% end

try % try a direct solver if memory is enough
    tmp_aaa=k\f;
catch
   parfor i=1:mesh.num_probes       
      tmp_aaa2(:,i)=pcg(k,f(:,i),1e-7,length(k),diag(diag(k))); 
   end
    
end
tmp_aaa=full(tmp_aaa).';

end









%/* ********************************** save data ********** */




function fem=save_data(itr,input,mesh,fem,aaa)
     

if itr<2 
       fem.geofac2=zeros(input.num_mes,1);
end

%/* SAVE ARRAY*/ 

for i=1:input.num_mes
	am=0; bm=0; bn=0; an=0;
	iim=mesh.pm(i); iin=mesh.pn(i);
	
    im=mesh.mes_nodes(iim);
	ia=mesh.pa(i); 
    
    if input.m_array_type==3
        in=0;
        ib=0;
    end
    
    if input.m_array_type==2
        in=mesh.mes_nodes(iin);
        ib=0;
    end
    
    if input.m_array_type==1
        in=mesh.mes_nodes(iin);
        ib=mesh.pb(i);
    end
    

	if((mesh.pa(i)*mesh.pm(i))~=0) ;am=aaa(ia,im); end
	if((mesh.pb(i)*mesh.pm(i))~=0) ;bm=aaa(ib,im); end
	if((mesh.pa(i)*mesh.pn(i))~=0) ;an=aaa(ia,in); end
	if((mesh.pb(i)*mesh.pn(i))~=0) ;bn=aaa(ib,in); end

	fem.array_model_data(i)=am-bm-an+bn;
% 	if (itr<2 )
%        fem.geofac2(i)=mesh.mean_res/fem.array_model_data(i) ;
% %       if (input.data_type==2)
% %          input.real_data(i)=input.real_data(i).*cm.geofac2(i);      
% %       end
%     end
    fem.array_model_data(i)=fem.array_model_data(i)*mesh.geofac(i);
%     fem.array_model_data(i)=fem.array_model_data(i)*fem.geofac2(i);


    
end

    for i=1:input.num_mes
        if fem.array_model_data(i)<0
               i 
               disp('error');
        elseif imag(fem.array_model_data(i))<0
              i 
               disp('quadrature error');
%                pause
        end
    end

%	/*  normalizing  data */
	if itr<2 
	  fem.normalizing_data=real(mesh.mean_res)./real(fem.array_model_data);
    end


 fem.array_model_data=fem.normalizing_data.*fem.array_model_data;

if itr==1 ;fem.array_model_data2=fem.array_model_data; end
end




