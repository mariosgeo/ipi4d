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


function [mesh]=initcm1(input,mesh,cx,cy,cd1,cd2)

 disp(' ')
 disp('======================')
 disp('=    ENTER INITCM1   =')
 disp('======================')
 disp(' ')

 %% tuning parameters
 factor = 2; % this is just for how many image samples you will use to calculate a tensor, the larger the factors, the less samples are used
 threshold_discontinuity=0.9  %0.987
 threshold_coherence=0.9934   %not used

 threshold_a1=22.5
 threshold_a2=67.5
 weight1=4      %3.85
 weight2=0.005  %0.05

 plof=2    %ra/factor
 %plof2=10

 p0=1.0
 p1=3.0   %1.0   %3.0

 plsf2_1=1
 plsf2_2=4
 %%end tuning parameters


 %% define size of training image
 % /!\ opposite conventions as SU ones
 n2=input.n1_TI
 n1=input.n2_TI


 %% define size of current model
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

 %% read image
 %% image=importdata(' ');
 fid=fopen(input.training_image,'r');
 %  int=fread(fid,[n2,n1],'float','b');
    image=fread(fid,[n2,n1],'float');
 fclose(fid);
 %  'b' means big endian. BTW java uses big endia


 %
 % INTERPOLATE MODEL ON IRREGULAR GRID
 %
 if input.interp_TI==1

    disp('INTERPOLATE TRAINING IMAGE ON IRREGULAR GRID')

    % input vectors
    vx1=[0:input.n2_TI-1]*input.hstep;    %row
    vz1=[0:input.n1_TI-1]'*input.hstep;   %column
    % output vectors
    vx2=mesh.XI(1,:);
    vz2=abs(mesh.YI(:,1));

    % complete initial vectors if output exceeds initial bounds
    if min(vx2)<min(vx1)
       vx1=[min(vx2) vx1];
       image=[image(:,1) image];
    end
    if max(vx2)>max(vx1)
       vx1=[vx1 max(vx2)];
       image=[image image(:,end)];
    end
    if min(vz2)<min(vz1)
       vz1=[min(vz2);vz1];
       image=[image(1,:);image];
    end
    if max(vz2)>max(vz1)
       vz1=[vz1;max(vz2)];
       image=[image;image(end,:)];
    end

    % interpolate
    image=interp2(vx1,vz1,image,vx2,vz2);
    size_image=size(image)
    size_map_param =size(mesh.map_param)

    % redefine sizes
    [n2,n1]=size(image)

 end   %end if interpolate training image


 %%
 %r1=n1/m1; r2=n2/m2;
 %ra=ceil(max(r1,r2));
 %plof=ra/factor


% import Java classes
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;

 % Constructs a filter with an isotropic Gaussian window
 lof = LocalOrientFilter(plof);

 %% Sets half-width of Gaussian derivative filter used to compute gradients
 %lof.setGradientSmoothing(4);

 %% Estimate 2-D structure tensors
 et = lof.applyForTensors(image);

 %% Inverts the structure tensors
 % p0 emphasizes overall amplitude and p1 emphasizes linearity (ie anisotropy).
 % For amplitude-independent tensors with all eigenvalues av equal to one, set p0 = 0.0.
 % To enhance linearity, set p1 > 1.0. To simply invert (and normalize) these tensors, set p0 = p1 = 1.0.
 et.invertStructure(p0,p1); % The are the p0, p1 in my paper
 % % the larger the second parameter, the more anisotropic the tensors,
 % % can be set to 3 or 4, but for geological cross-section image please leave as 1.0

 % ttm = aigi.Tensors(.1,.1,.1,image)
 % et = ttm.getTensors();

 %% Computes local semblance images using local smoothing filters. Local semblance (Hale, 2009) is defined to be a squared smoothed-image divided by a smoothed squared-image, where smoothing is performed by local smoothing filters along the eigenvectors of a structure tensor field.
 % Ref: Hale, D., 2009, Structure-oriented smoothing and semblance, CWP-635
 lsf1 = LocalSemblanceFilter(1,1);   %not used
 lsf2 = LocalSemblanceFilter(plsf2_1,plsf2_2);

 %% 2D smoothing directions correspond to eigenvectors of tensors.
 % The direction U corresponds to the largest eigenvalue (perpendicular to linear features),
 % The direction V corresponds to the smallest eigenvalue (parallel to linear features).
 V = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','V');
 U = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','U');

 %% Compute local semblance for the 2D array.
 % V,U = direction(s) for the first inner smoothing
 %  et = eigen-decomposition of the tensor field
 % image=array of input values
 s1 = lsf1.semblance(V,et,image);   %not used
 s2 = lsf2.semblance(U,et,image);


 %% IMAGE-GUIDED INTERPOLATION
 %??
 %s1=Sampling(n1,1.0,0);
 %s2=Sampling(n2,1.0,0);

 %?? syn_k not defined ??
 %bg=BlendedGridder2(et,syn_k,X,Y);
 %finalk=bg.grid(s1,s2);


 %% 
 X=ceil((mesh.param_x-min(mesh.param_x))*(n1/width));
 % Y=ceil((mesh.param_y2-min(mesh.param_y2))*(n2/depth));
 Y=ceil((mesh.param_y-min(mesh.param_y))*(n2/depth));

 %% with topography data, use mesh.param_y2
 X(find(X==0))=1;
 Y(find(Y==0))=1;
 X(find(X==n1))=n1-1;
 Y(find(Y==n2))=n2-1;
 X(find(X==n1+1))=n1-1;
 Y(find(Y==n2+1))=n2-1;
 YY=reshape(Y,m1,m2)';

 %% init. matrices
 w1=eye(m2*m1);
 w2=w1;
 w3=w1;
 w4=w1;

 Edge=zeros(n2,n1);
 Coherence=Edge;

 %% loop over columns
 for n=1:m1
     [A,B]=sort(s2(:,X(n)));
     [C,D]=sort(s1(:,X(n)));
     [m,J]=min(abs(A-threshold_discontinuity)); % This is the threshold that less than it is discontinuty
     Edge(:,X(n))=[B(1:J);zeros((n2-J),1)];
     [c,K]=min(abs(C-threshold_coherence)); % This is the threshold that larger than it is coherence
     Coherence(:,X(n))=[D(K:length(D));zeros((n2-length(D)+K-1),1)];
     J=unique(Edge(:,X(n)));

     K=unique(Coherence(:,X(n))); 
     for i=2:length(J)
         [x,j]=min(abs(J(i)-YY(:,n)));
         if (YY(j)>J(i))
            J(i)=j-1;
         else 
            J(i)=j;
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

     %   The same as above, if you have topography and want to avoid its influence
     %   if (K(i) == 1)%||K(i)==2)
     %      continue
     %   end

         i=(K(k)-1)*m1+n;
                    
         k1 = X(i); if k1==0; k1=1; elseif k1==n1-1 k1=n1-2; end
         k2 = Y(i); if k2==0; k2=1; elseif k2==n2-1 k2=n2-2; end

         % Gets the eigenvector v for the tensor at specified indices.
         % (returns the array {v1,v2} of eigenvector components)
         v1 = et.getEigenvectorV(k1,k2);
         v2 = et.getEigenvectorV(k1,k2+1);
         v3 = et.getEigenvectorV(k1+1,k2);
         v4 = et.getEigenvectorV(k1,k2-1);
         v5 = et.getEigenvectorV(k1-1,k2);
         a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/5)/pi*180;

         w1(i,i)=weight1;
         if (a<=threshold_a1 && a>=-threshold_a1)
            %w1(i,i)=40;
            w1(i,i)=weight1;
         elseif (a>threshold_a2 || a<-threshold_a2)
            if (n~=1 && n~=m1 && all(K(k)~=(m2-1:m2)) && all(K(k)~=(1:2)))
               w2(i,i)=weight1;
%           else 
%              w1(i,i)=4;
            end
                  
         elseif (a<-threshold_a1 && a>=-threshold_a2)
            if (all(K(k)~=(m2-1:m2)) && all(K(k)~=(1:2)))
               w3(i,i)=weight1;
%           else 
%              w1(i,i)=40;
            end

         elseif (a<=threshold_a2 && a>threshold_a1)
            if (all(K(k)~=(m2-1:m2)) && all(K(k)~=(1:2))) 
               w4(i,i)=weight1;
%           else 
%              w1(i,i)=40;
            end
          
         end   %end if threshold

     end   %end for k=2:length(K)


     for k=2:length(J)
%         % here if have topography and want to avoid treating layers close to
%         % air-ground boundary as discontinuty, uncommon the code below, can 
%         % choose to avoid first 1 or 2 or 3 or 4 or ... layers, this example avoids the first three layers 
%           if ((J(k) == 1))%%||(J(i)==2 )||(J(i)==3))%||(J(i)==4))
%             continue
%           end

         if n==1 
            n=2; 
         elseif n==m1 
            n=m1-1; 
         end

         i=(J(k))*m1+n; 
         k1 = X(i);if k1==0; k1=1; elseif k1==n1-1 k1=n1-2; end
         k2 = Y(i);if k2==0; k2=1; elseif k2==n2-1 k2=n2-2;  end

         v1 = et.getEigenvectorV(k1,k2);
         v2 = et.getEigenvectorV(k1,k2+1);
         v3 = et.getEigenvectorV(k1+1,k2);
         v4 = et.getEigenvectorV(k1,k2-1);
         v5 = et.getEigenvectorV(k1-1,k2);
         a = atan((v1(2)/v1(1)+v2(2)/v2(1)+v3(2)/v3(1)+v4(2)/v4(1)+v5(2)/v5(1))/5)/pi*180;
%        a = atan((v1(2)/v1(1)))/pi*180;

         if (a<=threshold_a1 && a>=-threshold_a1)
%           w1(i,i)=10;
            w2(i,i)=weight2;
            w3(i,i)=weight2;
            w4(i,i)=weight2;

         elseif (a<-threshold_a2 || a>threshold_a2)
%           w2(i,i)=10; 
            if J(k)~=m2 && n~=1 && n~=m1
               w1(i,i)=weight2;
               w3(i,i)=weight2;
               if i>m1
                  w4(i-(m1-1),i-(m1-1))=weight2; 
               end
            end

         elseif (a<-threshold_a1 && a>=-threshold_a2)  
            if J(k)~=m2 && n~=1
               w1(i,i)=weight2;
%              w3(i,i)=10;
               if i>m1
                  w2(i-m1,i-m1)=weight2;
               end
               if i>m1 
                  w4(i-(m1-1),i-(m1-1))=weight2;
               end
            end

         elseif (a<=threshold_a2 && a>threshold_a1)  
%           w4(i,i)=10;
            if J(k)~=m2 &&  n~=m1
               w1(i,i)=weight2;
               w2(i,i)=weight2;
               if i>(m1+1)
                  w3(i-(m1+1),i-(m1+1))=weight2;
               end 
            end
         end   %end if threshold
     end   %end for k=2:length(J)


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
%             if (a<=threshold_a1 && a>=-threshold_a1)
%                     w1(i,i)=200;
%             elseif (a>threshold_a2 || a<-threshold_a2)
%                     w2(i,i)=20;
%             elseif (a<-threshold_a1 && a>=-threshold_a2)
%                     w3(i,i)=20; 
%             elseif (a<=threshold_a2 && a>threshold_a1)
%                     w4(i,i)=20;
%           
%             end
%     end

 end   %end for n=1:m1


 % build covariance matrix
 ctc1=cx'*w1'*w1*cx;  
 ctc2=cy'*w2'*w2*cy;
 ctc3=cd1'*w3'*w3*cd1;
 ctc4=cd2'*w4'*w4*cd2;
 mesh.ctc = ctc1+ctc2+ctc3+ctc4; % this is the four-direction smoothing matrix

 % plot ellipses representing tensors
 ellipse

end   %end function

