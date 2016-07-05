%
% FUNCTION ACB: compute a locally-varying Lagrangian multiplier distribution  based on 
%               the data resolution, and applies it to the model covariance matrix.
%
% Author: Marios Karaoulis, Colorado School of Mines

% Modified by Francois Lavoue', Colorado School of Mines, September 2015.
% Modifications include:
% 1) some explicit comments for more readability,
% 2) some change of notations, in particular
%    - 'final' should be denoted 'inv' because it gathers all inversion variables
% 3) some modifications in
%    - the computation of the GN Hessian J'J (now including the data weighting matrix)
%    - the computation of the resolution matrix, which reflects only the data resolution

function [ctc,L1]=acb(input,mesh,fem)

 %FL: include data weighting matrix in the Gauss-Newton Hessian
 %JTJ=fem.array_jacobian.'*fem.array_jacobian;
 JTJ=fem.array_jacobian'*input.Wd*fem.array_jacobian;

 %FL: compute DATA resolution matrix only
 %    (the aim is to weight the model term according to data resolution)
 %resolution_matrix=(JTJ + 0.005*mesh.ctc)\(fem.array_jacobian.'*fem.array_jacobian);
 resolution_matrix=(JTJ + 0.005*eye(mesh.num_param))\JTJ;

 
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
 %(apply local weights from Lagrangian distribution)
 for i=1:mesh.num_param
     for j=1:mesh.num_param
         ctc(i,j)=mesh.ctc(i,j)*real(L1(i)); %real in cases of complex numbers
     end
 end



%--------------------PLOT LAGRANGIAN DISTRIBUTION------------------------
 if input.acb_plot==1
    %interpolate
    ZI=griddata(mesh.param_x,-mesh.param_y,real(L1),mesh.XI,mesh.YI);
    ZII = interp2(mesh.XI,mesh.YI,ZI,mesh.xxx,mesh.yyy);

    %plot
    figure;
    acb_fig1=gca;
    pcolor(mesh.xxx,mesh.yyy,ZII);

    %tune figure
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
 end

% Keep in for later plots

end   %end function acb

