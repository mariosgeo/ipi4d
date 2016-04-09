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
 

function [mesh,ctc,s1,s2]=initcm1(mesh,input,cx,cy,cd1,cd2)
%% 

 %% tuning parameters
 factor = 2; % this is just for how many image samples you will use to calculate a tensor, the larger the factors, the less samples are used
 threshold_discontinuity=0.9  %0.987
 threshold_coherence=0.9934   %not used

 threshold_a1=22.5
 threshold_a2=67.5
 weight1=20    %3.85
 weight2=10    %0.05

 plof=2    %ra/factor
 %plof2=10

 p0=1.0
 p1=3.0   %1.0   %3.0

 pslf2_1=1    %1
 pslf2_2=4 
 %%end tuning parameters


 %define sizes
 % /!\ opposite conventions as SU ones
 n2=input.n1_TI
 n1=input.n2_TI

 %depth=max(mesh.depth_n)
 %width=max(mesh.add_x_points)
 depth=max(mesh.param_y)
 width=max(mesh.param_x)

 %% image=importdata(' ');
 fid=fopen(input.training_image,'r');
 %  int=fread(fid,[n2,n1],'float','b');
    image=fread(fid,[n2,n1],'float');
 fclose(fid);
 %  'b' means big endian. BTW java uses big endian.

 %re-define sizes
 %WARNING: n1, n2 is the size of the guiding image,
 %         m1, m2 is the size of the current model
 %[n2,n1]=size(image);
 %mesh.m2=input.n2;
 %mesh.m1=input.n1;
 m2=mesh.m2
 m1=mesh.m1


 %% size ratio between guiding image and current model
 r1=n1/m1
 r2=n2/m2
 ra=ceil(max(r1,r2))

 %% import JTK classes
import interp.*
import edu.mines.jtk.dsp.*;

 %%
 lof = LocalOrientFilter(plof);
 et = lof.applyForTensors(image);
 et.invertStructure(p0,p1); % These are the p0, p1 in my paper
 % the larger the second parameter, the more anisotropic the tensors,
 % can be set to 3 or 4, but for geological cross-section image please leave as 1.0
 %

 lsf1 = LocalSemblanceFilter(16,4);   %not used
 lsf2 = LocalSemblanceFilter(pslf2_1,pslf2_2);   
 V = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','V');
 U = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','U');
 s1 = lsf1.semblance(V,et,image);
 s2 = lsf2.semblance(U,et,image);
 %

 %% convert indices in training image to indices in current model
 X=ceil(mesh.param_x*(n1/width));
 Y=ceil(mesh.param_y*(n2/depth)); %% with topography data, use mesh.param_y2
 X(find(X==0))=1;
 Y(find(Y==0))=1;
 YY=reshape(Y,m1,m2)';
 %

 %%
 w1=sparse(eye(m2*m1));
 w2=w1;
 w3=w1;
 w4=w1;

 Edge=zeros(n2,n1);
 Coherence=Edge;


 disp('ENTER IN LOOP')
 for n=1:m1
    [A,B]=sort(s2(:,X(n)));
    [C,D]=sort(s1(:,X(n)));  %not used
    [m,J]=min(abs(A-threshold_discontinuity)); % This is the threshold that less than it is discontinuty
    Edge(:,X(n))=[B(1:J);zeros((n2-J),1)];
    [c,K]=min(abs(C-threshold_coherence)); % This is the threshold that larger than it is coherence
    Coherence(:,X(n))=[D(K:length(D));zeros((n2-length(D)+K-1),1)];
    J=unique(Edge(:,X(n)));
    K=unique(Coherence(:,X(n))); 

    for i=2:length(J)
        [x,j]=min(abs(J(i)-YY(:,n)));
        J(i)=j;
    end

    for i=2:length(K)
        [x,K(i)]=min(abs(K(i)-YY(:,n)));
    end

    J=unique(J);
    K=unique(K);
    for i=2:length(J)-1
%         % here if have topography and want to avoid treating layers close to
%         % air-ground boundary as discontinuty, uncommon the code below, can 
%         % choose to avoid first 1 or 2 or 3 or 4 or ... layers, this example avoids the first three layers 
%           if ((J(i) == 1)||(J(i)==2 )||(J(i)==3))%||(J(i)==4))
%             continue
%           end
        i=(J(i)-1)*m1+n;
        k1 = X(i);
        k2 = Y(i);

        %compute average of eigen vectors
        val_sum=0;
        if k1<m1 && k2<m2
           v1 = et.getEigenvectorV(k1,k2);
           val_sum=val_sum+1;
        else
           v1=[1 ; 0];
        end

        if k1<m1 && k2+1<m2
           v2 = et.getEigenvectorV(k1,k2+1);
           val_sum=val_sum+1;
        else
           v2=[1 ; 0];
        end

        if k1+1<m1 && k2<m2
           v3 = et.getEigenvectorV(k1+1,k2);
           val_sum=val_sum+1;
        else
           v3=[1 ; 0];
        end  

        if k1<m1
           v4 = et.getEigenvectorV(k1,k2-1);
           val_sum=val_sum+1;
        else
           v4=[1 ; 0];
        end

        if k2<m2
           v5 = et.getEigenvectorV(k1-1,k2);
           val_sum=val_sum+1;
        else
           v5=[1 ; 0];
        end

        %compute angle
        a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/val_sum)/pi*180;

        if a<=threshold_a1 && a>=-threshold_a1   %if a<=22.5 && a>=-22.5
                    w1(i,i)=weight1;             %   w1=200   %3.85
                    w2(i,i)=weight2;             %   w2=0.05
                    w3(i,i)=weight2;
                    w4(i,i)=weight2;

        elseif a<-threshold_a2 || a>threshold_a2   %if a<-67.5 || a>67.5
                    w2(i,i)=weight1;
                    w1(i,i)=weight2;
                    w3(i,i)=weight2;
                    if i>m1
                        w4(i-(m1-1),i-(m1-1))=weight2; 
                    end

        elseif a<-threshold_a1 && a>=-threshold_a2  
                    w1(i,i)=weight2;
                    w3(i,i)=weight1;
                    if i>m1
                        w2(i-m1,i-m1)=weight2;
                    end
                    if i>m1 
                        w4(i-(m1-1),i-(m1-1))=weight2;
                    end

        elseif a<=threshold_a2 && a>threshold_a1  
                    w4(i,i)=weight1;
                    w1(i,i)=weight2;
                    w2(i,i)=weight2;
                    if i>(m1+1)
                        w3(i-(m1+1),i-(m1+1))=weight2;
                    end         
        end   %end if threshold

    end   %end for i=2:length(J)-1


%      % These are code for coherence but I usually don't use them...
% for i=2:length(K)
% % %    The same as above, if you have topography and want to avoid its influence
% %     if (K(i) == 1)%||K(i)==2)
% %             continue
% %       end
%         i=(K(i)-1)*m1+n;
%         k1 = X(i);
%         k2 = Y(i);
%         v1 = et.getEigenvectorV(k1,k2);
%         v2 = et.getEigenvectorV(k1,k2+1);
%         v3 = et.getEigenvectorV(k1+1,k2);
%         v4 = et.getEigenvectorV(k1,k2-1);
%         v5 = et.getEigenvectorV(k1-1,k2);
%         a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/5)/pi*180;   
%             if a<=22.5 && a>=-22.5
%                     w1(i,i)=8;
%             elseif a>67.5 || a<-67.5
%                     w2(i,i)=8;
%             elseif a<-22.5 && a>=-67.5
%                     w3(i,i)=8; 
%             else a<=67.5 && a>22.5
%                     w4(i,i)=8;
%           
%             end
%     end

        
 end   %end for n=1:m1


 % compute covariance matrix
 ctc1=cx'*w1'*w1*cx;  %plot the diag to see if it contains discontinuities in the y-dir
 ctc2=cy'*w2'*w2*cy;  %  "        "                                               x-dir
 ctc3=cd1'*w3'*w3*cd1;
 ctc4=cd2'*w4'*w4*cd2;
 ctc = ctc1+ctc2+ctc3+ctc4; % this is the four-direction smoothing matrix


 % plot tensors
 ellipse


end   %end function


