function [ctc,L1]=acb(input,mesh,fem)

global acb_fig1 


   JTJ=fem.array_jacobian.'*fem.array_jacobian;

   %ctc=c'*c;
   resolution_matrix=(JTJ + 0.005*mesh.ctc)\(fem.array_jacobian.'*fem.array_jacobian);



 
 
% -----------------Calculate SP matrix-------------------------------------
   for i=1:mesh.num_param
       cur1=i;
       tmp999=0.0;
       for j=1:mesh.num_param
           cur2=j;
           if (cur1>mesh.num_param) ;cur1=i-mesh.num_param; end
           if (cur2>mesh.num_param) ;cur2=j-mesh.num_param; end
           
           	x1_dis=mesh.tmp_param(cur1,1);
    		x2_dis=mesh.tmp_param(cur2,1);
		
        	y1_dis=mesh.tmp_param(cur1,2);
            y2_dis=mesh.tmp_param(cur2,2);
	    	weight=( (x1_dis - x2_dis)*(x1_dis - x2_dis) + (y1_dis - y2_dis)*(y1_dis - y2_dis) ); %/*Not sqrt because in SP i need weiht*weight*/
           	
            tmp999=tmp999+ weight*( (1-mesh.S(i,j))*resolution_matrix(i,j) ) * ( (1-mesh.S(i,j))*resolution_matrix(i,j) ); 
       end
       
       SP(i)=tmp999;
   end

SP_max=max(SP);
SP_min=min(SP);


% -----------Security check in case of low resolution----------------------
% L1
minimum_security=0.5;
SP_security=minimum_security*( tmp999/mesh.num_param);

for i=1:mesh.num_param
    tmp999=0;
    for j=1:mesh.num_param
        tmp999=tmp999+resolution_matrix(i,j);
    end
    if (tmp999 <minimum_security) SP(i)=SP_max; end
end







% ----------------Lagrange Distribution------------------------------------
for i=1:mesh.num_param
	L1(i)=log10(input.lagrn_min) +  ( (log10(input.lagrn_max) - log10(input.lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(SP(i)) - log10(SP_min));    
    L1(i)=10^(L1(i));
end  



ctc=zeros(mesh.num_param,mesh.num_param);
%update now ctc for export
for i=1:mesh.num_param
    for j=1:mesh.num_param
        ctc(i,j)=mesh.ctc(i,j)*real(L1(i)); %real in cases of complex numbers
    end
end













%--------------------PLOT LAGRANGIANG DISTRIBUTION------------------------

    ZI=griddata(mesh.param_x,-mesh.param_y,real(L1),mesh.XI,mesh.YI);
    ZII = interp2(mesh.XI,mesh.YI,ZI,mesh.xxx,mesh.yyy);
    pcolor(acb_fig1,mesh.xxx,mesh.yyy,ZII);
    colorbar('peer',acb_fig1);
    shading(acb_fig1, 'interp');
    title(acb_fig1,'Lagrangian Distribution');
    xlabel(acb_fig1,'Distance (m)');
    ylabel(acb_fig1,'Depth (m)');
    axis(acb_fig1,'equal');
    xlim(acb_fig1,[mesh.min_x mesh.max_x]);
    ylim(acb_fig1,[-mesh.max_y -mesh.min_y]);
    %colorbar(acb_fig1);
    drawnow;

% Keep in for later plots







end