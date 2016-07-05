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
% instead of doing 
%   dX = (J'J + beta*Cm)^(-1) * J'*Cd*[G(X)-d)]
% do this:
%   dX = reshape(interp.Sm.getLaplacianI(mesh.et2,reshape(J'*Cd*[G(X)-d)],mesh.m1,mesh.m2)',(J'J)',scale,mesh.s)',mesh.num_param,1);
% Here scale is how you want to scale the tensors, the larger the more emphasize on structure
% mesh.s is using structur-oriented semblance s1 to scale the tensors, if don't wanna use can just put [] instead of mesh.s
% 


function [mesh,s2]=initcm2(mesh,input)

 %% tuning parameters
 factor = 2; % this is just for how many image samples you will use to calculate a tensor, the larger the factors, the less samples are used

 scale=1
 mesh.scale=scale;

 plof=2   %ra/factor

 p0=1.0
 p1=3.0   %1.0   %3.0

 plsf1_1=16   %not used
 plsf1_2=4

 plsf2_1=1    %4
 plsf2_2=4    %16
 %%end tuning parameters


 %% define size of training image
 % /!\ opposite conventions as SU ones
 n2=input.n1_TI
 n1=input.n2_TI


 %% define size of current model
 [m2,m1]=size(mesh.map_param);
 mesh.m2=m2;
 mesh.m1=m1;
 mesh.num_param=m1*m2

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



%
% JAVA LIBRARY
%
import interp.*
import edu.mines.jtk.dsp.*;


 % Constructs a filter with an isotropic Gaussian window
 lof = LocalOrientFilter(plof);

 %% Sets half-width of Gaussian derivative filter used to compute gradients
 % lof.setGradientSmoothing(2);

 %% Estimate 2-D structure tensors
 et = lof.applyForTensors(image);

 %% Inverts the structure tensors
 % p0 emphasizes overall amplitude and p1 emphasizes linearity (ie anisotropy).
 % For amplitude-independent tensors with all eigenvalues av equal to one, set p0 = 0.0.
 % To enhance linearity, set p1 > 1.0. To simply invert (and normalize) these tensors, set p0 = p1 = 1.0.
 et.invertStructure(p0,p1); % The are the p0, p1 in my paper
 % the larger the second parameter, the more anisotropic the tensors,
 % can be set to 3 or 4, but for geological cross-section image please leave as 1.0

 %% Computes local semblance images using local smoothing filters. Local semblance (Hale, 2009) is defined to be a squared smoothed-image divided by a smoothed squared-image, where smoothing is performed by local smoothing filters along the eigenvectors of a structure tensor field.
 % Ref: Hale, D., 2009, Structure-oriented smoothing and semblance, CWP-635
 lsf1 = LocalSemblanceFilter(plsf1_1,plsf1_2);   %not used
 lsf2 = LocalSemblanceFilter(plsf2_1,plsf2_2);

 %% 2D smoothing directions correspond to eigenvectors of tensors.
 % The direction U corresponds to the largest eigenvalue (perpendicular to linear features),
 % The direction V corresponds to the smallest eigenvalue (parallel to linear features).
 V = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','V');
 U = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','U');
 s1 = lsf1.semblance(V,et,image);
 s2 = lsf2.semblance(U,et,image);
 mesh.et2 = EigenTensors2(m1,m2);

 %% 
 X=ceil((mesh.param_x-min(mesh.param_x))*(n1/width));
 Y=ceil((mesh.param_y-min(mesh.param_y))*(n2/depth));
 %% with topography data, use mesh.param_y2
 % Y=ceil((mesh.param_y2-min(mesh.param_y2))*(n2/depth));

 X(find(X==0))=1;
 Y(find(Y==0))=1;
 X(find(X==n1))=n1-1;
 Y(find(Y==n2))=n2-1;
 YY=reshape(Y,m1,m2)';

 %% build scaling matrix using semblance
 mesh.s = zeros(m2,m1);
 ev = interp.Sm.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));

 for i=1:mesh.num_param

     %define k1
     k1 = X(i);
     if X(i)==X(1)
        k1= X(2);
     else 
        if X(i)==X(m1)
           k1= X(m1-1);
        end
     end

     %define k2
     k2 = Y(i);
     if Y(i)==YY(1)
        k2= YY(2);
     else 
        if Y(i)==YY(m2)
           k2= YY(m2-2);
        end
     end

     %define i1, i2
     i1 = mod(i,m1);
     i2 = ceil(i/m1)-1;
     if (i1==0) 
        i1=m1-1; 
     else 
        i1=i1-1; 
     end

%     a = et.getTensor(k1,k2);
%     mesh.et2.setTensor(i1,i2,a);
%     mesh.et2.invertStructure(1,1);
%     a = et.getEigenvalues(k1,k2);
     a = ev(:,k2,k1);
     aa = et.getEigenvectorU(k1,k2);
%     if abs(a(1)-a(2))<0.99
%        a(1)=0.25; a(2)=1;
%         aa(1)=0;aa(2)=1;
%     end
     a(2)=1;

     mesh.et2.setEigenvalues(i1,i2,a);
     mesh.et2.setEigenvectorU(i1,i2,aa);
     mesh.s(i2+1,i1+1)=1/(1.0001-s2(k2,k1));
 end   %end for i=1:mesh.num_param

 %% plot ellipses
 ellipse

end   %end function
