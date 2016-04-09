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
% 
% There are some parameters that you can tune, please see the comments.
% 


function [mesh,ctc]=initcm1(mesh,cx,cy,cd1,cd2)

%% 
factor = 3; % this is just for how many image samples you will use to calculate a tensor, the larger the factors, the less samples are used

%% define sizes
[m2,m1]=size(mesh.map_param);
mesh.m2=m2;
mesh.m1=m1; 

depth=max(mesh.param_y);
width=max(mesh.param_x);
%% alternate definitions
% depth=max(mesh.depth_n); 
% depth = depth - min(mesh.param_y2);
% width=max(mesh.add_x_points);
% depth=mesh.max_y+mesh.min_y;
% depth = depth - min(mesh.param_y2);

% % m=mesh.X;m(7)=m(6);m(6)=33;
% m=[18,39,21,39,60,12,30];
% slope=(m(7)-m(6))/m(5);
% a=m(6)+m(1)*slope;
% b=m(6)+m(3)*slope;
% c=m(6)+m(2)*slope;
% d=m(6)+m(4)*slope;
% figure
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
% xlim([-9 width])
% % set(gcf,'units','normalized','outerposition',[0 0 1 1])
%  %%
% [XX,map]=getframe(gca);
% XX=XX(2:end-1,2:end-1,:);
% % imwrite(XX,'ch_im.pgm','pgm');
% clear XX


% image=importdata('fttest.pgm');
%  fid=fopen(' .dat','r');
%  int=fread(fid,[n2,n1],'float','b');
%  fclose(fid);
%  'b' means big endian. BTW java uses big endian.
image=importdata('ch_im.pgm');
[n2,n1]=size(image);
% load('amplitude.mat')
% maxf=max((amplitude(:,1)));
% % minf=min(log10(amplitude(:,1)));
% for i=1:4
%     f(i)=((amplitude(i,1)))/maxf*255;
% end
% image(find(image==150))=f(1);
% image(find(image==226))=f(2);
% image(find(image==29))=f(3);
% image(find(image==76))=f(4);
% image(find(image==0))=10;
 image(find(image==150))=100;
%%
r1=n1/m1; r2=n2/m2;
ra=ceil(max(r1,r2));

% import interp.*
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
% javaaddpath('/home/jzhou/work/yma/trunk/build/jar/yma.jar')
%%
lof = LocalOrientFilter(8);
% lof.setGradientSmoothing(4);
et = lof.applyForTensors(image);
et.invertStructure(1,1); % The are the p0, p1 in my paper
% % the larger the second parameter, the more anisotropic the tensors,
% % can be set to 3 or 4, but for geological cross-section image please leave as 1.0
% 

% ttm = aigi.Tensors(.1,.1,.1,image)
% et = ttm.getTensors();
lsf1 = LocalSemblanceFilter(1,1);
lsf2 = LocalSemblanceFilter(1,1);
V = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','V');
U = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','U');
s1 = lsf1.semblance(V,et,image);
s2 = lsf2.semblance(U,et,image);

s1=Sampling(n1,1.0,0);
s2=Sampling(n2,1.0,0);
bg=BlendedGridder2(et,syn_k,X,Y);
finalk=bg.grid(s1,s2);

%% 
X=ceil((mesh.param_x-min(mesh.param_x))*(n1/width));
% Y=ceil((mesh.param_y2-min(mesh.param_y2))*(n2/depth));
Y=ceil((mesh.param_y-min(mesh.param_y))*(n2/depth));

%% with topography data, use mesh.param_y2
X(find(X==0))=1;Y(find(Y==0))=1;
X(find(X==n1))=n1-1;Y(find(Y==n2))=n2-1;
X(find(X==n1+1))=n1-1;Y(find(Y==n2+1))=n2-1;
YY=reshape(Y,m1,m2)';

%%
w1=eye(m2*m1);w2=w1;w3=w1;w4=w1;
Edge=zeros(n2,n1);Coherence=Edge;
for n=1:m1
    [A,B]=sort(s2(:,X(n)));
    [C,D]=sort(s1(:,X(n)));
    [m,J]=min(abs(A-0.9)); % This is the threshold that less than it is discontinuty
    Edge(:,X(n))=[B(1:J);zeros((n2-J),1)];
    [c,K]=min(abs(C-0.9934)); % This is the threshold that larger than it is coherence
    Coherence(:,X(n))=[D(K:length(D));zeros((n2-length(D)+K-1),1)];
    J=unique(Edge(:,X(n)));
    
    K=unique(Coherence(:,X(n))); 
    for i=2:length(J)
        [x,j]=min(abs(J(i)-YY(:,n)));
        if (YY(j)>J(i))
            J(i)=j-1;
        else J(i)=j;
         end
      end
    for i=2:length(K)
        [x,K(i)]=min(abs(K(i)-YY(:,n)));
    end
    
    J=unique(J);
    K=unique(K);
  for l=2:length(J)
        if J(l)==J(l-1)+1;
            J(l)=J(l-1);
        end
    end
           J=unique(J); 
    for k=2:length(K)
% %    The same as above, if you have topography and want to avoid its influence
%     if (K(i) == 1)%||K(i)==2)
%             continue
%       end

        i=(K(k)-1)*m1+n;
                    
                    k1 = X(i);if k1==0; k1=1; elseif k1==n1-1 k1=n1-2; end
        k2 = Y(i);if k2==0; k2=1;elseif k2==n2-1 k2=n2-2;  end
        v1 = et.getEigenvectorV(k1,k2);
        v2 = et.getEigenvectorV(k1,k2+1);
        v3 = et.getEigenvectorV(k1+1,k2);
        v4 = et.getEigenvectorV(k1,k2-1);
        v5 = et.getEigenvectorV(k1-1,k2);
        a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/5)/pi*180;
  w1(i,i)=4;
        if (a<=22.5 && a>=-22.5)
%     w1(i,i)=40;
   w1(i,i)=4;
elseif (a>67.5 || a<-67.5)
    if (n~=1 && n~=m1 && all(K(k)~=(m2-1:m2)) && all(K(k)~=(1:2))) w2(i,i)=4;
%                   else w1(i,i)=4;
                     end
                  
            elseif (a<-22.5 && a>=-67.5)
                     if (all(K(k)~=(m2-1:m2)) && all(K(k)~=(1:2))) w3(i,i)=4;
%                     else w1(i,i)=40;
                     end
            elseif (a<=67.5 && a>22.5)
                    if (all(K(k)~=(m2-1:m2)) && all(K(k)~=(1:2))) w4(i,i)=4;
%                     else w1(i,i)=40;
                    end
          
            end

    end
     for k=2:length(J)
%         % here if have topography and want to avoid treating layers close to
%         % air-ground boundary as discontinuty, uncommon the code below, can 
%         % choose to avoid first 1 or 2 or 3 or 4 or ... layers, this example avoids the first three layers 
%           if ((J(k) == 1))%%||(J(i)==2 )||(J(i)==3))%||(J(i)==4))
%             continue
%           end
if n==1 n=2; elseif n==m1 n==m1-1; end
    i=(J(k))*m1+n; 
  k1 = X(i);if k1==0; k1=1; elseif k1==n1-1 k1=n1-2; end
        k2 = Y(i);if k2==0; k2=1;elseif k2==n2-1 k2=n2-2;  end
        v1 = et.getEigenvectorV(k1,k2);
        v2 = et.getEigenvectorV(k1,k2+1);
        v3 = et.getEigenvectorV(k1+1,k2);
        v4 = et.getEigenvectorV(k1,k2-1);
        v5 = et.getEigenvectorV(k1-1,k2);
        a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/5)/pi*180;
%         a = atan((v1(2)/v1(1)))/pi*180;

        if (a<=22.5 && a>=-22.5)
%                     w1(i,i)=10;
                    w2(i,i)=0.005;
                    w3(i,i)=0.005;
                    w4(i,i)=0.005;
            elseif (a<-67.5 || a>67.5)
%                     w2(i,i)=10; 
                    if J(k)~=m2 && n~=1 && n~=m1
                    w1(i,i)=0.005;
                    w3(i,i)=0.005;
                    if i>m1
                        w4(i-(m1-1),i-(m1-1))=0.005; 
                    end
                    end
            elseif (a<-22.5 && a>=-67.5)  
                    if J(k)~=m2 && n~=1
                w1(i,i)=0.005;
%                     w3(i,i)=10;
                    if i>m1
                        w2(i-m1,i-m1)=0.005;
                    end
                    if i>m1 
                        w4(i-(m1-1),i-(m1-1))=0.005;
                    end
                    end
            elseif (a<=67.5 && a>22.5)  
%                     w4(i,i)=10;
                    if J(k)~=m2 &&  n~=m1
                    w1(i,i)=0.005;
                    w2(i,i)=0.005;
                    if i>(m1+1)
                        w3(i-(m1+1),i-(m1+1))=0.005;
                    end 
                    end
        end
     end
    
     
     % These are code for coherence but I usually don't use them...
% for k=2:length(K)
% % %    The same as above, if you have topography and want to avoid its influence
% %     if (K(i) == 1)%||K(i)==2)
% %             continue
% %       end
%         i=(K(k)-1)*m1+n;
%         k1 = X(i);
%         k2 = Y(i);
%         v1 = et.getEigenvectorV(k1,k2);
% %         v2 = et.getEigenvectorV(k1,k2+1);
% %         v3 = et.getEigenvectorV(k1+1,k2);
% %         v4 = et.getEigenvectorV(k1,k2-1);
% %         v5 = et.getEigenvectorV(k1-1,k2);
% %         a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/5)/pi*180;
%    a = atan((v1(2)/v1(1)))/pi*180;
%             if (a<=22.5 && a>=-22.5)
%                     w1(i,i)=200;
%             elseif (a>67.5 || a<-67.5)
%                     w2(i,i)=20;
%             elseif (a<-22.5 && a>=-67.5)
%                     w3(i,i)=20; 
%             elseif (a<=67.5 && a>22.5)
%                     w4(i,i)=20;
%           
%             end
%     end

        
end
  ctc1=cx'*w1'*w1*cx;  
  ctc2=cy'*w2'*w2*cy;
  ctc3=cd1'*w3'*w3*cd1;
  ctc4=cd2'*w4'*w4*cd2;
  ctc = ctc1+ctc2+ctc3+ctc4; % this is the four-direction smoothing matrix
ellipse
end
