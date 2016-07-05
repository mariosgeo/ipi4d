% Author: Jieyi Zhou, Colorado School of Mines, 2014
% mesh is just a structure carrying parameters around
% please pass your own parameter and replace them
% m2 = number of rows of model matrix
% m1 = number of columns of model matrix
% depth = total depth of model section
% width = total width of model section
% mesh.param_x is the x-axis location vector of each cell of model, row after row
% mesh.param_y is the y-axis location (depth) vector of each cell of model, row after row
% mesh.param_y2 is the depth vector of each cell of model with topography, row after row

% mesh.num_param is the number of total model cells, = m1*m2

% There are some parameters that you can tune, please see the comments.
% 
 
% How to use this is after initializing (calling initcm2 before inversion)
% instead of doing dX = (J'J + beta*Cm)^(-1) * J'*Cd*[G(X)-d)], do this:
% dX = reshape(interp.Sm.getLaplacianI(mesh.et2,reshape(J'*Cd*[G(X)-d)],mesh.m1,mesh.m2)',(J'J)',scale,mesh.s)',mesh.num_param,1);
% Here scale is how you want to scale the tensors, the larger the more emphasize on structure
% mesh.s is using structur-oriented semblance s1 to scale the tensors, if don't wanna use can just put [] instead of mesh.s
% 

function [mesh]=initcm2(mesh,scale)
%% 
factor = 2; % this is just for how many image samples you will use to calculate a tensor, the larger the factors, the less samples are used

[m2,m1]=size(mesh.map_param);
mesh.m2=m2; mesh.m1=m1; 
depth=mesh.max_y+mesh.min_y;%(mesh.depth_n); 
% depth = depth - min(mesh.param_y2);

width=42;%max(mesh.add_x_points);
mesh.scale=scale;

%%
% 
% % width = 127.5;
% % depth = 34.43;
% % % m=[4.965;11.57;6.962;14.23;20.37;33;57];
% % %   m=[4;10;6;15;19;33;50];
% % m(7)=m(6);m(6)=33;
% % slope=(m(7)-m(6))/m(5);
% % a=m(6)+m(1)*slope;
% % b=m(6)+m(3)*slope;
% % c=m(6)+m(2)*slope;
% % d=m(6)+m(4)*slope;
% figure
% patch([0,0,width,width],[m(1),depth,depth,m(1)],'b','LineStyle','None')
% patch([0;0;width;width;(m(5)+m(6))/2],[m(2),m(1),m(1),m(4),m(3)],'g','LineStyle','None')
% patch([0,0,(m(5)+m(6))/2,width,width],[0,m(2),m(3),m(4),0],'k','LineStyle','None')
% patch([30,m(5),m(6),55],[0,m(1),m(1),0],'r','LineStyle','None')
% % patch([0,0,a,m(6)],[0,m(1),m(1),0],'g','LineStyle','None')
% % patch([b;m(6);width;width],[m(3),0,0,m(3)],'b','LineStyle','None')
% % patch([0,0,c,a],[m(1),m(2),m(2),m(1)],'k','LineStyle','None')
% % patch([b,d,width,width],[m(3),m(4),m(4),m(3)],'g','LineStyle','None')
% % patch([0,0,width,width,m(7),c],[m(2),depth,depth,m(5),m(5),m(2)],'r','LineStyle','None')
% % patch([d,m(7),width,width],[m(4),m(5),m(5),m(4)],'k','LineStyle','None')

%  slope=(m(7)-m(6))/m(5);
% a=m(6)+m(1)*slope;
% b=m(6)+m(3)*slope;
% c=m(6)+m(2)*slope;
% d=m(6)+m(4)*slope;
% figure
% 
% patch([-9,-9,a,m(6)],[0,m(1),m(1),0],'y','LineStyle','None')
% patch([b;m(6);width;width],[m(3),0,0,m(3)],'g','LineStyle','None')
% patch([-9,-9,c,a],[m(1),m(2),m(2),m(1)],'b','LineStyle','None')
% patch([b,d,width,width],[m(3),m(4),m(4),m(3)],'y','LineStyle','None')
% patch([-9,-9,width,width,m(7),c],[m(2),depth,depth,m(5),m(5),m(2)],'r','LineStyle','None')
% patch([d,m(7),width,width],[m(4),m(5),m(5),m(4)],'b','LineStyle','None')
% 
% axis ij
% axis equal
% axis off
% ylim([0 depth])
% xlim([0 width])
% set(gcf,'units','normalized','outerposition',[0 0 1 1])
 %%
% [XX,map]=getframe(gca);
% XX=XX(2:end-1,2:end-1,:);
% imwrite(XX,'ch_im.pgm','pgm');
% clear XX
% %%
% image=importdata('p3test.pgm');
% image=importdata('seismic-2014pagosafinal.pgm');
image=importdata('ch_im.pgm');
% image=importdata('PAGO08.pgm');
[n2,n1]=size(image);
% image=image(2:n2-1,2:n1-1);
% [n2,n1]=size(image);
image(find(image==0))=10;
% close(gcf)

%%
r1=n1/m1; r2=n2/m2;
ra=ceil(max(r1,r2));

import interp.*
import edu.mines.jtk.dsp.*;


% import edu.mines.jtk.dsp.*;
% javaaddpath('/home/jzhou/work/yma/trunk/build/jar/yma.jar')
% 
% 
% ttm = aigi.Tensors(.1,.1,.1,image)
% et = ttm.getTensors();
%%
lof = LocalOrientFilter(2.5);
% lof.setGradientSmoothing(2);
et = lof.applyForTensors(image);
et.invertStructure(0.01,1); % The are the p0, p1 in my paper
% the larger the second parameter, the more anisotropic the tensors,
% can be set to 3 or 4, but for geological cross-section image please leave as 1.0
%%
lsf1 = LocalSemblanceFilter(4,1);
lsf2 = LocalSemblanceFilter(1,4);
V = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','V');
U = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','U');
s1 = lsf1.semblance(V,et,image);
s2 = lsf2.semblance(U,et,image);
mesh.et2 = EigenTensors2(m1,m2);

%% 
X=ceil((mesh.param_x-min(mesh.param_x))*(n1/width));
Y=ceil((mesh.param_y-min(mesh.param_y))*(n2/depth)); %% with topography data, use mesh.param_y2
% Y=ceil((mesh.param_y2-min(mesh.param_y2))*(n2/depth));
X(find(X==0))=1;Y(find(Y==0))=1;
X(find(X==n1))=n1-1;Y(find(Y==n2))=n2-1;
YY=reshape(Y,m1,m2)';
%%
mesh.s = zeros(m2,m1);
ev = interp.Sm.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));
for i=1:mesh.num_param
    k1 = X(i);
    if X(i)==X(1)
        k1= X(2);
    else if X(i)==X(m1)
            k1= X(m1-1);
        end
    end
    k2 = Y(i);
    if Y(i)==YY(1)
        k2= YY(2);
    else 
if Y(i)==YY(m2)
            k2= YY(m2-2);
        end
    end
    i1 = mod(i,m1);if(i1==0) i1=m1-1; else i1=i1-1; end
    i2 = ceil(i/m1)-1; 
%     a = et.getTensor(k1,k2);
%     mesh.et2.setTensor(i1,i2,a);
%     mesh.et2.invertStructure(1,1);
%     a = et.getEigenvalues(k1,k2);
a =  ev(:,k2,k1);
   
    aa = et.getEigenvectorU(k1,k2);
%     if abs(a(1)-a(2))<0.99
%        a(1)=0.25; a(2)=1;
%         aa(1)=0;aa(2)=1;
%     end
a(2)=1;
    mesh.et2.setEigenvalues(i1,i2,a);
    mesh.et2.setEigenvectorU(i1,i2,aa);
    mesh.s(i2+1,i1+1)=1/(1.0001-s2(k2,k1));
end
ellipse

end