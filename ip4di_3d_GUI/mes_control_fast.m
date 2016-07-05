function [fem,ex]=mes_control_fast(itr,input,mesh,fem,current_file)



   disp('--------- FEM Calculations --------');
   fem=mes_control4(itr,input,mesh,fem);
   
   [ex,fem,mesh]=save_data(itr,input,mesh,fem);

 
   disp('-------- JACOBIAN CALCULATION -------');
   
if input.time_lapse_flag==1;input.jacobian_flag=1;end  
   if(input.jacobian_flag==1 || itr==1) 
%         fem=jac_flux_fast2(itr,fem,mesh,input) ; % non parallel
        fem=jac_flux_fast2_par(itr,fem,mesh,input) ; % parallel
   end

   if(input.jacobian_flag==2 && itr>1) ;fem=quasi_newton(input,mesh,fem); end
   if(input.jacobian_flag==3 && (itr==4 || itr==7 || itr==10))
%         fem=jac_flux_fast2(itr,fem,mesh,input) ; % non parallel
        fem=jac_flux_fast2_par(itr,fem,mesh,input) ; % parallel 
   elseif (input.jacobian_flag==3 && itr~=1)
        fem=quasi_newton(input,mesh,fem);
   end






end





function fem=mes_control4(itr,input,mesh,fem)


fprintf('\n');
fprintf('%s %d %s ', 'Solving for ', mesh.num_probes, ' electrodes');
fprintf('\n');


% Here take analytical solution for homogeneous earth
if itr==1 
    fem=secondary_field(fem,mesh,input);
end

% Return the LU of the finite element matrix. Return the F matrix
[k,f,fem]=form(itr,mesh,input,fem);


if mesh.num_nodes<=20000
  fem.pot=k\f;% LU factorization, or Gaussian elimination 
  % fem.pot=fem.Ap'+k\((fem.Kp-k)*fem.Ap');  
else
   parfor i=1:mesh.num_probes
      pot(:,i)=pcg(k,f(:,i),1e-3,length(k),diag(diag(k)));
   end
fem.pot=pot;
end


fem.pot=full(fem.pot.');
end



function [k,f,fem]=form(itr,mesh,input,fem)


% issym=@(x) all(all(x==x.')); % Macro to test if symmetric
current=1;
%Keep kernel matrix
if itr==1
        fem.kernel=zeros(8,8,mesh.num_elements);
end

xc=0; yc=0; zc=0;


b_node=zeros(8,8);
c_node=zeros(8,8);
d_node=zeros(8,8);

% I am using k_index so I can avoid preallocate k matrix


a_index=zeros(8*8*mesh.num_elements,1);
b_index=zeros(8*8*mesh.num_elements,1);
val_index=zeros(8*8*mesh.num_elements,1);
ke=zeros(8,8);  



k_index=1;

for elno=1:mesh.num_elements

    if itr==1 % Calculate element dimesnions only ONCE == KERNEL FUNCTION         
    % find element's dimensions */

% %                 for i=1:8
% %                   for j=i:8
% %                          ii=mesh.icon(i,elno); jj=mesh.icon(j,elno);
% %                          a_2=abs(mesh.x_node_coord(ii)-mesh.x_node_coord(jj));
% %                          if(a_2~=0) ;xc=(mesh.x_node_coord(ii)+mesh.x_node_coord(jj))/2; break; end
% %                   end
% %                          if(a_2~=0) ;break; end
% %                 end
            ii=mesh.icon(:,elno);
            xs=unique(mesh.x_node_coord(ii));
            a_2=abs(xs(1)-xs(2));
            xc=(xs(1)+xs(2))/2;

            ys=unique(mesh.y_node_coord(ii));
            b_2=abs(ys(1)-ys(2));
            yc=(ys(1)+ys(2))/2;

            zs=unique(mesh.z_node_coord(ii));
            c_2=abs(zs(1)-zs(2));
            zc=(zs(1)+zs(2))/2;

%                for i=1:8
%                   for j=i:8
%                       ii=mesh.icon(i,elno); jj=mesh.icon(j,elno);
%                       b_2=abs(mesh.y_node_coord(ii)-mesh.y_node_coord(jj));
%                          if(b_2~=0) ;yc=(mesh.y_node_coord(ii)+mesh.y_node_coord(jj))/2; break; end
%                   end
%                          if(b_2~=0) ;break; end
%                end

%                for i=1:8
%                   for j=i:8
%                       ii=mesh.icon(i,elno); jj=mesh.icon(j,elno);
%                          c_2=abs(mesh.z_node_coord(ii)-mesh.z_node_coord(jj));
%                          if(c_2~=0) ;zc=(mesh.z_node_coord(ii)+mesh.z_node_coord(jj))/2; break; end
%                   end
%                          if(c_2~=0) ;break; end
%                end


                a=a_2/2; b=b_2/2; c=c_2/2;

            %/* set_up variables */
% iii=mesh.icon(:,elno);


xii=mesh.x_node_coord(ii)-xc;
yii=mesh.y_node_coord(ii)-yc;
zii=mesh.z_node_coord(ii)-zc;

xixj=xii*xii';
yiyj=yii*yii';
zizj=zii*zii';

b_node=((xixj)/(8*a*a*a*b*c)).* ...
                  (c*c*b*b+(c*c*(yiyj))/3+(b*b*(zizj))/3+((yiyj).*(zizj))/9);
              
c_node=((yiyj)/(8*a*b*b*b*c)).* ...
                  (c*c*a*a+(a*a*(zizj))/3+(c*c*(xixj))/3+((xixj).*(zizj))/9);  
              
              
d_node=((zizj)/(8*c*c*c*b*a)).* ...
                  (a*a*b*b+(b*b*(xixj))/3+(a*a*(yiyj))/3+((yiyj).*(xixj))/9);              

% %             for i=1:8
% %                for j=i:8
% %                   ii=mesh.icon(i,elno); jj=mesh.icon(j,elno);
% %                   xi=mesh.x_node_coord(ii)-xc;
% %                   xj=mesh.x_node_coord(jj)-xc;
% %                   yi=mesh.y_node_coord(ii)-yc;
% %                   yj=mesh.y_node_coord(jj)-yc;
% %                   zi=mesh.z_node_coord(ii)-zc;
% %                   zj=mesh.z_node_coord(jj)-zc;
% % 
% %                   b_node(i,j)=((xi*xj)/(8*a*a*a*b*c))* ...
% %                   (c*c*b*b+(c*c*yi*yj)/3+(b*b*zi*zj)/3+(yi*yj*zi*zj)/9);
% % 
% %                   c_node(i,j)=((yi*yj)/(8*a*b*b*b*c))* ...
% %                   (c*c*a*a+(a*a*zi*zj)/3+(c*c*xi*xj)/3+(xi*xj*zi*zj)/9);
% % 
% %                   d_node(i,j)=((zi*zj)/(8*c*c*c*b*a))* ...
% %                   (a*a*b*b+(b*b*xi*xj)/3+(a*a*yi*yj)/3+(yi*yj*xi*xj)/9);
% %                end
% %             end

%/* form element matrices*/

% ke=(b_node+c_node+d_node)*2./mesh.prop(elno);
ke=(b_node+c_node+d_node);
% for i=1:8      
%     for j=i:8
%         ke(i,j)=(fem.b_node(i,j)+fem.c_node(i,j)+fem.d_node(i,j))*2/mesh.prop(elno);
%     end
% end
% for i=1:8
%     for j=i:8
% 			ke(j,i)=ke(i,j);
%     end
% end

% ke = ke + ke.' - diag(diag(ke));
% Calculate element dimesnions only ONCE == KERNEL FUNCTION         
fem.kernel(:,:,elno)=ke;
end % End of calculation of element dimensions


ke=fem.kernel(:,:,elno)*2./mesh.prop(elno);




    
%assemble all element stiffness terms in to the global stiffness matrix */
for i=1:8
    ii=mesh.icon(i,elno);
        for j=1:8
            jj=mesh.icon(j,elno);
%           fem.k(ii,jj)=fem.k(ii,jj)+ke(i,j);
            a_index(k_index)=ii;
            b_index(k_index)=jj;
            val_index(k_index)=ke(i,j);
            k_index=k_index+1;   
            
        end
end



end % end of form




% Now create the k matrix
k=sparse(a_index,b_index,val_index,mesh.num_nodes,mesh.num_nodes);



    


% Impose the Dirichlet boundary conditions    
if itr==1
    fem.prop_ebc=zeros(mesh.num_ebc,mesh.num_probes);
%      for j=1:mesh.num_probes
%            tmp11=sqrt((mesh.x_node_coord(mesh.node_ebc)-mesh.probe_x(j)).*(mesh.x_node_coord(mesh.node_ebc)-mesh.probe_x(j))+ ...
%                        (mesh.y_node_coord(mesh.node_ebc)-mesh.probe_y(j)).*(mesh.y_node_coord(mesh.node_ebc)-mesh.probe_y(j))+ ...
%                        (mesh.z_node_coord(mesh.node_ebc)-mesh.probe_z(j)).*(mesh.z_node_coord(mesh.node_ebc)-mesh.probe_z(j)));
%     	   
%             tmp12=sqrt((mesh.x_node_coord(mesh.node_ebc)-mesh.probe_x(j)).*(mesh.x_node_coord(mesh.node_ebc)-mesh.probe_x(j))+ ...
%     			       (mesh.y_node_coord(mesh.node_ebc)-mesh.probe_y(j)).*(mesh.y_node_coord(mesh.node_ebc)-mesh.probe_y(j))+ ...
%     			       (-mesh.z_node_coord(mesh.node_ebc)-mesh.probe_z(j)).*(-mesh.z_node_coord(mesh.node_ebc)-mesh.probe_z(j)));
%             fem.prop_ebc(:,j)=mesh.mean_res./(tmp11.*tmp12.*4.*pi./(tmp11+tmp12));
%      end
    
    

%     tmp11=sqrt((mesh.x_node_coord(mesh.node_ebc)-mesh.xc_mesh).*(mesh.x_node_coord(mesh.node_ebc)-mesh.xc_mesh)+ ...
%             (mesh.y_node_coord(mesh.node_ebc)-mesh.yc_mesh).*(mesh.y_node_coord(mesh.node_ebc)-mesh.yc_mesh)+ ...
% 		     mesh.z_node_coord(mesh.node_ebc).*mesh.z_node_coord(mesh.node_ebc));
%     fem.prop_ebc=mesh.mean_res./(tmp11.*2.*pi);
%     
%     fem.prop_ebc=zeros(mesh.num_ebc,1);     
     
end








f=zeros(mesh.num_nodes,mesh.num_probes);



%--------------------First Attempt ---------------------------------------
% f2=zeros(mesh.num_nodes,mesh.num_probes);
% k2=full(k2);
% for t=1:mesh.num_probes
%     t
%     k3=k2;
%     for i=1:mesh.num_ebc
%        f2(mesh.node_ebc(i),t)=fem.prop_ebc(i,t);
%         k3(mesh.node_ebc(i),mesh.node_ebc(i))=1;
%         for j=1:mesh.num_nodes
% 		      if (mesh.node_ebc(i)==j); continue; end
%                   f2(j,t)=f2(j,t)-k3(j,mesh.node_ebc(i))*fem.prop_ebc(i,t);
%                   k3(mesh.node_ebc(i),j)=0;
%                   k3(j,mesh.node_ebc(i))=0;
%               
%         end
%     end    
% end

%---------------------Second Attempt --------------------------------------





%---------------Accurate part----------------------------------------------


% % for t=1:mesh.num_probes
% %     for i=1:mesh.num_ebc
% % %        f(mesh.node_ebc(i),t)=prop_ebc(i,t); 
% % % %        jlist = setdiff(1:mesh.num_nodes,mesh.node_ebc(i));  
% %        rem=mesh.node_ebc(i);
% %        jlist=[1:rem-1 rem+1:mesh.num_nodes;];
% %        f(jlist,t)= f(jlist,t)-sum(k2(jlist,mesh.node_ebc(i)).*fem.prop_ebc(i,t),2);
% %        
% %     end 
% % end


%-------------------Fast part----------------------------------------------

num_nodes=mesh.num_nodes;
prop_ebc=fem.prop_ebc;
node_ebc=mesh.node_ebc;

for t=1:mesh.num_probes  % do nout use parfor here. running out of memory
%     prop2=ones(num_nodes,1)*prop_ebc(:,t).';
%     f(:,t)=f(:,t)-sum(k(:,node_ebc).*prop2,2);   
%     f(:,t)=f(:,t)-sum(k(:,node_ebc)*diag(prop_ebc(:,t)),2);   
%     f(:,t)=f(:,t)-sum(  k(:,node_ebc).*repmat(prop_ebc(:,t),size(k(:,node_ebc),1),1),2);   
%      f(:,t)=f(:,t)-sum( bsxfun(@times,k(:,node_ebc),prop_ebc(:,t).') ,2);  
     

    f(:,t)=f(:,t)-sum( bsxfun(@times,k(:,node_ebc),prop_ebc(:,t).') ,2);

end

k(node_ebc,:)=0;
k(:,node_ebc)=0;
k(sub2ind(size(k),mesh.node_ebc,mesh.node_ebc))=1;  
f(mesh.node_ebc,:)=fem.prop_ebc;


% % Second attempt
% % k(sub2ind(size(k),mesh.node_ebc,mesh.node_ebc))=10e70;  
% % f(mesh.node_ebc,:)=10e70*fem.prop_ebc;


% % create source's here
for t=1:mesh.num_probes
    source=mesh.mes_nodes(t);
    f(source,t)=f(source,t)+2*current;
%     f(mesh.sink,t)=f(mesh.cinc,t)-2*current;
end

% keep stiffness matrix from homogeneous earth 
if itr==1  fem.Kp=k;
            for t=1:mesh.num_probes
                fem.Ap(t,mesh.node_ebc)=fem.prop_ebc(:,t);
            end
end 



f=sparse(f);
k=sparse(k);

% fem.k=k;
% fem.f=f;
end



function fem=secondary_field(fem,mesh,input)

% Here I calculate the potentials in case of homogeneous earth using
% analytical solution and store it in matrix Ap
current=2;
fem.Ap=zeros(mesh.num_probes,mesh.num_nodes);

val=real(mesh.mean_res);

for i=1:mesh.num_probes
        
            tmp11=sqrt((mesh.x_node_coord-mesh.probe_x(i)).*(mesh.x_node_coord-mesh.probe_x(i))+ ...
                       (mesh.y_node_coord-mesh.probe_y(i)).*(mesh.y_node_coord-mesh.probe_y(i))+ ...
                       (mesh.z_node_coord-mesh.probe_z(i)).*(mesh.z_node_coord-mesh.probe_z(i)));
            tmp12=sqrt((mesh.x_node_coord-mesh.probe_x(i)).*(mesh.x_node_coord-mesh.probe_x(i))+ ...
                       (mesh.y_node_coord-mesh.probe_y(i)).*(mesh.y_node_coord-mesh.probe_y(i))+ ...
                       (-mesh.z_node_coord-mesh.probe_z(i)).*(-mesh.z_node_coord-mesh.probe_z(i)));
            fem.Ap(i,:)=val./(tmp11.*tmp12.*4.*pi./(tmp11+tmp12));
       
end
% singularity removal
for i=1:mesh.num_probes
    j=mesh.mes_nodes(i);  
        tmp13=0.1;
        tmp11=sqrt((tmp13)*(tmp13)+(tmp13)*(tmp13)+(tmp13)*(tmp13));
        tmp12=sqrt((tmp13)*(tmp13)+(tmp13)*(tmp13)+ ...
                (-mesh.z_node_coord(j)-mesh.probe_z(i))*(-mesh.z_node_coord(j)-mesh.probe_z(i)));
        fem.Ap(i,j)=val./(tmp11*tmp12*4*pi/(tmp11+tmp12));
end



% for i=1:mesh.num_probes
%     j=mesh.mes_nodes(i);
%     nums=[1:j-1 j+1:mesh.num_nodes];
%     tmp=fem.Ap(nums);
%     xs=mesh.x_node_coord(nums);
%     ys=mesh.y_node_coord(nums);
%     
%     tmp1=TriScatteredInterp(xs,ys,real(tmp)');
%     tmp1=tmp1(mesh.x_node_coord(j),mesh.y_node_coord(j));
%     
%     tmp2=TriScatteredInterp(xs,ys,imag(tmp)');
%     tmp2=tmp2(mesh.x_node_coord(j),mesh.y_node_coord(j));
%     
%     fem.Ap(i,j)=complex(tmp1,tmp2);
%     
% 
% end



end



function [ex,fem,mesh]=save_data(itr,input,mesh,fem)

ex=0;
if itr==1
    fem.pot_norm=zeros(mesh.num_probes,mesh.num_nodes);
    fem.array_model_data=zeros(input.num_mes,1);
    %fem.pot_model_data=zeros(input.num_mes,1);
    fem.oldemes=zeros(input.num_mes,1);
    fem.norm_data=zeros(input.num_mes,1);
    mesh.geofac2=zeros(input.num_mes,1);
end
%normalizing  data*/
if itr==1 
    fem.pot_norm=fem.Ap./real(fem.pot);
%     [ind,ind2]=find(fem.pot_norm==Inf);
%     fem.pot_norm(isnan(fem.pot_norm))=0;
%     fem.pot_norm(sub2ind(size(fem.pot_norm),ind,ind2))=0;
%      fem.pot_norm(isinf(fem.pot_norm)|| isnan(fem.pot_norm))=0;
  fem.pot_norm(~isfinite(fem.pot_norm))=0;
    
    fem.array_model_data2=fem.array_model_data;
end


% % %perform normalization */
fem.pot=fem.pot.*fem.pot_norm;


% Array model data
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
	
	if((mesh.pa(i)*mesh.pm(i))~=0) am=fem.pot(ia,im); end
	if((mesh.pb(i)*mesh.pm(i))~=0) bm=fem.pot(ib,im); end 
	if((mesh.pa(i)*mesh.pn(i))~=0) an=fem.pot(ia,in); end
	if((mesh.pb(i)*mesh.pn(i))~=0) bn=fem.pot(ib,in); end
	
	fem.array_model_data(i)=am-bm-an+bn;
    if itr==1
       mesh.geofac2(i)=mesh.mean_res/fem.array_model_data(i) ;
      if (input.data_type==2)
         input.real_data(i)=input.real_data(i).*mesh.geofac2(i);      
      end
    end
%     fem.array_model_data(i)=fem.array_model_data(i)*mesh.geofac2(i);
 	fem.array_model_data(i)=fem.array_model_data(i)*mesh.geofac(i);%
    
end

% if( (itr==1 && input.bgr_res_flag==0) || (itr==0 && input.bgr_res_flag~=0) )
%     for i=1:input.num_mes ;fem.norm_data(i)=mesh.mean_res/fem.array_model_data(i);end
% end    
%      fem.array_model_data=fem.array_model_data.*fem.norm_data;
 

if itr==1 ;fem.oldmes=fem.array_model_data; end



end

function fem=jac_flux_fast2(itr,fem,mesh,input)

fem.array_jacobian=zeros(input.num_mes,mesh.num_param);
% The main idea is that for all elements that belonh to a paramter create a
% vectrors aa,ab,am,an where they have the potentials 

[ii,jj]=meshgrid(1:mesh.num_param,1:input.num_mes);

for idx = 1:numel(ii) 
    
    ia=mesh.pa(jj(idx)); im=mesh.pm(jj(idx)); ib=mesh.pb(jj(idx)); in=mesh.pn(jj(idx));
 
    
    jtmp_this_iter =nonzeros( mesh.jtmp(:,ii(idx))); 
    
        
    c=fem.kernel(:,:,jtmp_this_iter); 

    tmp=reshape(mesh.icon(:,jtmp_this_iter),1,8*length(jtmp_this_iter)); 

    
    N=length(jtmp_this_iter);

    if(mesh.pa(jj(idx))==0)
        c_ka=zeros(8,N);
    else
        ka=reshape(fem.pot(ia,tmp),8,[],N);
        c_ka=squeeze(  sum(   bsxfun(@times,c,ka) ,   1)  );
    end
    
    if(mesh.pb(jj(idx))==0)
        c_kb=zeros(8,N);
    else
        kb=reshape(fem.pot(ib,tmp),8,[],N);
        c_kb=squeeze(  sum(   bsxfun(@times,c,kb) ,   1)  ); 
    end
    
    if(mesh.pa(jj(idx))==0)
        km=zeros(8,N);
    else
        km=reshape(fem.pot(im,tmp),8,N);
    end
    
    if(mesh.pa(jj(idx))==0)
        kn=zeros(8,N);
    else
        kn=reshape(fem.pot(in,tmp),8,N);
    end
    
      
    all=sum( (c_ka-c_kb).*(km-kn));        
    sum2=sum(all)/(mesh.prop(jtmp_this_iter(1))^2);
 
%     fem.array_jacobian(jj(idx),ii(idx))=sum2*mesh.geofac2(jj(idx));
    fem.array_jacobian(jj(idx),ii(idx))=sum2*mesh.geofac(jj(idx));
    %/*find log jac */
    fem.array_jacobian(jj(idx),ii(idx))=fem.array_jacobian(jj(idx),ii(idx))*mesh.res_param1(ii(idx))/(fem.array_model_data(jj(idx))*input.stdev_error(jj(idx)));
    
end 





end



function fem=jac_flux_fast2_par(itr,fem,mesh,input)

fem.array_jacobian=zeros(input.num_mes,mesh.num_param);

% The main idea is that for all elements that belonh to a paramter create a
% vectrors aa,ab,am,an where they have the potentials 

[ii,jj]=meshgrid(1:mesh.num_param,1:input.num_mes);
array_jacobian=zeros(numel(ii),1);


% create here extra variable to avoid overhead
jtmp=mesh.jtmp;
kernel=fem.kernel;
pot=fem.pot;
array_model_data=fem.array_model_data;
stdev_error=input.stdev_error;
res_param1=mesh.res_param1;
pa=mesh.pa;
pb=mesh.pb;
pm=mesh.pm;
pn=mesh.pn;
icon=mesh.icon;
prop=mesh.prop;
geofac=mesh.geofac;

parfor idx = 1:numel(ii) 
    
    ia=pa(jj(idx)); im=pm(jj(idx)); ib=pb(jj(idx)); in=pn(jj(idx));
 
    jtmp_this_iter =nonzeros( jtmp(:,ii(idx))); 
        
    c=kernel(:,:,jtmp_this_iter); 
    tmp=reshape(icon(:,jtmp_this_iter),1,8*length(jtmp_this_iter)); 
    
    N=length(jtmp_this_iter);

    if(pa(jj(idx))==0)
        c_ka=zeros(8,N);
    else
        ka=reshape(pot(ia,tmp),8,[],N);
        c_ka=squeeze(  sum(   bsxfun(@times,c,ka) ,   1)  );
    end
    
    if(pb(jj(idx))==0)
        c_kb=zeros(8,N);
    else
        kb=reshape(pot(ib,tmp),8,[],N);
        c_kb=squeeze(  sum(   bsxfun(@times,c,kb) ,   1)  ); 
    end
    
    if(pm(jj(idx))==0)
        km=zeros(8,N);
    else
        km=reshape(pot(im,tmp),8,N);
    end
    
    if(pn(jj(idx))==0)
        kn=zeros(8,N);
    else
        kn=reshape(pot(in,tmp),8,N);
    end
    
      
    all=sum( (c_ka-c_kb).*(km-kn));        
    sum2=sum(all)/(prop(jtmp_this_iter(1))^2);
 

    array_jacobian(idx)=sum2;
    %/*find log jac */  
    
end 


%now create array_jacobian
for idx = 1:numel(ii) 
   fem.array_jacobian(jj(idx),ii(idx))=array_jacobian(idx)*geofac(jj(idx))*res_param1(ii(idx))/(array_model_data(jj(idx))*stdev_error(jj(idx)));    
end


end













function fem=jac_flux(itr,fem,mesh,input)
% global fem cm input
% b_node=zeros(8,8);
% c_node=zeros(8,8);
% d_node=zeros(8,8);

tic
fem.array_jacobian=zeros(input.num_mes,mesh.num_param);
% fem.array_jacobian2=zeros(input.num_mes,mesh.num_param);
jtmp=zeros(1500,1);



%dlg = ProgressDialog();


% kernel=fem.kernel;
% pot=fem.pot;
% icon=mesh.icon;
% pa=mesh.pa;pb=mesh.pb;pm=mesh.pm;pn=mesh.pn;





for i=1:mesh.num_param
%   disp(i)
    %pause(0.1)
    % update progress bar
   % dlg.FractionComplete = i/mesh.num_param;
    % update status message
   % dlg.StatusMessage = sprintf('%d%% complete', fix(100*i/mesh.num_param)); 
    
    m=0;
    for k=1:mesh.num_elements
        if(mesh.elem_param(k)==i) ;m=m+1; jtmp(m)=k;end
    end
    
    
    for j=1:input.num_mes
        jam=0;    jbm=0; jan=0; jbn=0;
        ia=mesh.pa(j); ib=mesh.pb(j); im=mesh.pm(j); in=mesh.pn(j);
        sum2=0;
                
        for n=1:m
            el=jtmp(n);
            
% % % % % % %                 for i1=1:8
% % % % % % %                   for j1=i:8
% % % % % % %                          ii=mesh.icon(i1,el); jj=mesh.icon(j1,el);
% % % % % % %                          a_2=abs(mesh.x_node_coord(ii)-mesh.x_node_coord(jj));
% % % % % % %                          if(a_2~=0) ;xc=(mesh.x_node_coord(ii)+mesh.x_node_coord(jj))/2; break; end
% % % % % % %                   end
% % % % % % %                          if(a_2~=0) ;break; end
% % % % % % %                 end
% % % % % % % 
% % % % % % %                for i1=1:8
% % % % % % %                   for j1=i:8
% % % % % % %                       ii=mesh.icon(i1,el); jj=mesh.icon(j1,el);
% % % % % % %                       b_2=abs(mesh.y_node_coord(ii)-mesh.y_node_coord(jj));
% % % % % % %                          if(b_2~=0) ;yc=(mesh.y_node_coord(ii)+mesh.y_node_coord(jj))/2; break; end
% % % % % % %                   end
% % % % % % %                          if(b_2~=0) ;break; end
% % % % % % %                end
% % % % % % % 
% % % % % % %                for i1=1:8
% % % % % % %                   for j1=i:8
% % % % % % %                       ii=mesh.icon(i1,el); jj=mesh.icon(j1,el);
% % % % % % %                          c_2=abs(mesh.z_node_coord(ii)-mesh.z_node_coord(jj));
% % % % % % %                          if(c_2~=0) ;zc=(mesh.z_node_coord(ii)+mesh.z_node_coord(jj))/2; break; end
% % % % % % %                   end
% % % % % % %                          if(c_2~=0) ;break; end
% % % % % % %                end
% % % % % % % 
% % % % % % % 
% % % % % % %                 a=a_2/2; b=b_2/2; c=c_2/2;
% % % % % % % 
% % % % % % %             %/* set_up variables */
% % % % % % % 
% % % % % % %             for i1=1:8
% % % % % % %                for j1=1:8
% % % % % % %                   ii=mesh.icon(i1,el); jj=mesh.icon(j1,el);
% % % % % % %                   xi=mesh.x_node_coord(ii)-xc;
% % % % % % %                   xj=mesh.x_node_coord(jj)-xc;
% % % % % % %                   yi=mesh.y_node_coord(ii)-yc;
% % % % % % %                   yj=mesh.y_node_coord(jj)-yc;
% % % % % % %                   zi=mesh.z_node_coord(ii)-zc;
% % % % % % %                   zj=mesh.z_node_coord(jj)-zc;
% % % % % % % 
% % % % % % %                   b_node(i1,j1)=((xi*xj)/(8*a*a*a*b*c))* ...
% % % % % % %                   (c*c*b*b+(c*c*yi*yj)/3+(b*b*zi*zj)/3+(yi*yj*zi*zj)/9);
% % % % % % % 
% % % % % % %                   c_node(i1,j1)=((yi*yj)/(8*a*b*b*b*c))* ...
% % % % % % %                   (c*c*a*a+(a*a*zi*zj)/3+(c*c*xi*xj)/3+(xi*xj*zi*zj)/9);
% % % % % % % 
% % % % % % %                   d_node(i1,j1)=((zi*zj)/(8*c*c*c*b*a))* ...
% % % % % % %                   (a*a*b*b+(b*b*xi*xj)/3+(a*a*yi*yj)/3+(yi*yj*xi*xj)/9);
% % % % % % %                end
% % % % % % %             end
% % % % % % %             
% % % % % % %             c=b_node+c_node+d_node;
            c=fem.kernel(:,:,el);

            tmp=mesh.icon(1:8,el);
            
%             ka=fem.pot(ia,tmp);
%             kb=fem.pot(ib,tmp);
%             km=fem.pot(im,tmp);
%             kn=fem.pot(in,tmp);
            
            
            
            if((mesh.pa(j)*mesh.pm(j))~=0)
                %jam=calc_flx(el,ia,im);
                tmp=mesh.icon(1:8,el);
                a1(1:8)=fem.pot(ia,tmp);%Keep in mind that a1'=conj(a1)
                a2(1:8)=fem.pot(im,tmp);
                S2=((conj(a1))'*a2).*c;
%                 S2=(a1.'*a2).*c;
                
%                 S2=(ka.'*km).*c;                
                jam=sum(S2(:));
            end
            if((mesh.pb(j)*mesh.pm(j))~=0) 
                %jbm=calc_flx(el,ib,im);
                %tmp=mesh.icon(1:8,el);
                a1(1:8)=fem.pot(ib,tmp);%Keep in mind that a1'=conj(a1)
                a2(1:8)=fem.pot(im,tmp);
                S2=((conj(a1))'*a2).*c;
%                 S2=(a1.'*a2).*c;
                
%                 S2=(kb.'*km).*c;  
                jbm=sum(S2(:));
          end
            if((mesh.pa(j)*mesh.pn(j))~=0) 
                %jan=calc_flx(el,ia,in);%Keep in mind that a1'=conj(a1)
                tmp=mesh.icon(1:8,el);
                a1(1:8)=fem.pot(ia,tmp);
                a2(1:8)=fem.pot(in,tmp);
                S2=((conj(a1))'*a2).*c;
%                 S2=(a1.'*a2).*c;
                
%                 S2=(ka.'*kn).*c;  
                jan=sum(S2(:));
            end
            if((mesh.pb(j)*mesh.pn(j))~=0) 
                %jbn=calc_flx(el,ib,in);%Keep in mind that a1'=conj(a1)
                tmp=mesh.icon(1:8,el);
                a1(1:8)=fem.pot(ib,tmp);
                a2(1:8)=fem.pot(in,tmp);
                S2=((conj(a1))'*a2).*c;
%                 S2=(a1.'*a2).*c;
                
%                 S2=(kb.'*kn).*c;  
                jbn=sum(S2(:));
            end
            sum2=sum2+(jam-jbm-jan+jbn)./(mesh.prop(el)^2);
        end
%     fem.array_jacobian(j,i)=sum2*mesh.geofac(j);%*2/mesh.probe_spacing_x;
  fem.array_jacobian(j,i)=sum2*mesh.geofac(j);

    %/*find log jac */
    fem.array_jacobian(j,i)=fem.array_jacobian(j,i)*mesh.res_param1(i)/fem.array_model_data(j);
    fem.array_jacobian(j,i)=fem.array_jacobian(j,i)/input.stdev_error(j);
    end
end
%delete(dlg);
toc
end

function sum2=calc_flx(elm,cr,pt)

global cm fem
el=elm;  sum2=0;
m1=cr;
m2=pt;

ke=zeros(2,8);
%/* find element's dimensions */        
   for i=1:8
      for j=i:8
             ii=mesh.icon(i,el); jj=mesh.icon(j,el);
			 a_2=abs(mesh.x_node_coord(ii)-mesh.x_node_coord(jj));
			 if(a_2~=0)  ;xc=(mesh.x_node_coord(ii)+mesh.x_node_coord(jj))/2; break;end
      end
			 if(a_2~=0) break; end
   end
                 
   for i=1:8
      for j=i:8
             ii=mesh.icon(i,el); jj=mesh.icon(j,el);
			 b_2=abs(mesh.y_node_coord(ii)-mesh.y_node_coord(jj));
			 if(b_2~=0)  yc=(mesh.y_node_coord(ii)+mesh.y_node_coord(jj))/2; break; end 
	  end	 
             if(b_2~=0) ;break;end
      
   end
             
   for i=1:8
      for j=i:8
             ii=mesh.icon(i,el); jj=mesh.icon(j,el);
			 c_2=abs(mesh.z_node_coord(ii)-mesh.z_node_coord(jj));
			 if(c_2~=0)  zc=(mesh.z_node_coord(ii)+mesh.z_node_coord(jj))/2; break;end
      end		 
             if(c_2~=0) ;break;end        
         end

a=a_2/2; b=b_2/2; c=c_2/2;
%/* set_up variables */   
   
for i=1:8
 
  for j=1:8
     
	 ii=mesh.icon(i,el); jj=mesh.icon(j,el);
     xi=mesh.x_node_coord(ii)-xc;
     xj=mesh.x_node_coord(jj)-xc;
     yi=mesh.y_node_coord(ii)-yc;
     yj=mesh.y_node_coord(jj)-yc;
     zi=mesh.z_node_coord(ii)-zc;
     zj=mesh.z_node_coord(jj)-zc;
     
	 b_node(i,j)=((xi*xj)/(8*a*a*a*b*c))* ...
     (c*c*b*b+(c*c*yi*yj)/3+(b*b*zi*zj)/3+(yi*yj*zi*zj)/9);
     
     c_node(i,j)=((yi*yj)/(8*a*b*b*b*c))* ...
     (c*c*a*a+(a*a*zi*zj)/3+(c*c*xi*xj)/3+(xi*xj*zi*zj)/9);

     d_node(i,j)=((zi*zj)/(8*c*c*c*b*a))* ...
     (a*a*b*b+(b*b*xi*xj)/3+(a*a*yi*yj)/3+(yi*yj*xi*xj)/9);

  end
end

for k=1:8
   tmp=mesh.icon(k,el);
   ke(1,k)=fem.pot(m1,tmp);
end

if(m1~=m2)
   for k=1:8
	 tmp=mesh.icon(k,el);
	 ke(2,k)=fem.pot(m2,tmp);
   end

 for i=1:8
   for j=1:8	
	 sum2=sum2+(ke(1,i)*ke(2,j)*b_node(i,j)+ ...
	 	   ke(1,i)*ke(2,j)*c_node(i,j)+ ...
		   ke(1,i)*ke(2,j)*d_node(i,j)) ...
		  /(mesh.prop(el)*mesh.prop(el));
   end
 end
end


end