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
%
% There are some parameters that you can tune, please see the comments.
%
 
% Modified by F. Lavoue', Colorado School of Mines, September 2015
% in order to:
% 1) consider TI and grid of different sizes.
% 2) comment the use of Dave Hale's and Jieyi Zhou's Java routines
%    (comments are based on Dave Hale's documentation at http://dhale.github.io/jtk/api)

%JZ: How to use this is after initializing (calling initcm2 before inversion)
% instead of doing in invert_cntr.m
%   dX = (J'J + beta*Cm)^(-1) * J'*Cd*[G(X)-d)]
% do this:
%   dX = reshape(interp.Sm.getLaplacianI(mesh.et2,reshape(J'*Cd*[G(X)-d)],mesh.m1,mesh.m2)',(J'J)',mesh.scale,mesh.s,input.lagrn)',mesh.num_param,1);
% Here mesh.scale is how you want to scale the tensors, the larger the more emphasize on structure.
%   FL: this is redundant with the role of the Lagrangian multiplier input.lagrn, better don't use mesh.scale.
%   FL2JZ: can you confirm my understanding?
%   Answer: YES
% mesh.s is using structure-oriented semblance s2 to (locally) scale the tensors, if don't wanna use can just put [] instead of mesh.s
%   FL: this is redundant with the spatially-varying Lagrangian defined by ACB (but slightly different, though...)
%   FL2JZ: can you confirm my understanding?
%   Answer: YES, if we define Cm^-1 = G'DG, but not exactly if Cm^-1 = I + G'DG...


function [mesh,s2]=initcm2(input,mesh)

disp(' ')
disp('========================================================')
disp('=                   ENTER INITCM2                      =')
disp('= define structure-based directional Laplacian filters =')
disp('========================================================')
disp(' ')

 %%=== TUNING PARAMETERS ===%%
 % FL: should be defined in 12_inversion_parameters.m
 %     to avoid hard-coding here...

 %factor = 2;   %JZ: this is just for how many image samples you will use to calculate a tensor,
 %                   the larger the factors, the less samples are used
 %FL: not used anymore

 mesh.scale=1;   %FL: should probably not be used anymore (ie let to 1),
                 %    because it is redundant with the spatially varying Lagrangian distribution defined by ACB.

 plof=4;   % half-width of isotropic Gaussian window for LocalOrientFilter 
           % FL2JZ: is it expressed in nb of cells in TI?
           % Answer: YES (should be at least 4, but may depend on the nb of TI cells in 1 inversion cell)
           % NB:it is also the diffusivity alpha used for interpolation (see DH's CWP report)

 %JZ: These are the p0, p1 in my paper.
 %    The larger the second parameter, the more anisotropic the tensors,
 %    can be set to 3 or 4, but for geological cross-section image please leave as 1.0.
 p0=0.0;
 p1=1.0;   %3.0
 %DH's doc for invertStructure:
 % p0 emphasizes overall amplitude and p1 emphasizes linearity (ie anisotropy).
 % For amplitude-independent tensors with all eigenvalues av equal to one, set p0 = 0.0.
 % To enhance linearity, set p1 > 1.0. To simply invert (and normalize) these tensors, set p0 = p1 = 1.0.


 %% Parameters of local smoothing filters (lsf) for local semblance computation
 %plsf1_1=16;   %not used
 %plsf1_2=4;

 plsf2_1=1;    %4    % half-width of 1st smoothing filter
 plsf2_2=4;    %16   % half-width of 2nd smoothing filter
 %%=== END TUNING PARAMETERS ===%%


 %% define size of training image
 % /!\ opposite conventions as SU ones
 n2=input.n1_TI;
 n1=input.n2_TI;
 vx_TI=[0:input.hstep:(n1-1)*input.hstep];
 vy_TI=[0:input.hstep:(n2-1)*input.hstep];

 %% read image
 fid=fopen(input.training_image,'r');
    image=fread(fid,[n2,n1],'float');
 fclose(fid);

 %% define size of current model
 [m2,m1]=size(mesh.map_param);
 mesh.m1=m1;
 mesh.m2=m2;

 if input.debug==1
 % define domain size
 % (just to make sure it is the same on both grids, or not necessarily...)
    depth=max(vy_TI)
    width=max(vx_TI)
    depth=max(mesh.param_y)
    width=max(mesh.param_x)
 end

 %% alternate definitions
 % depth=max(mesh.depth_n); 
 % depth = depth - min(mesh.param_y2);
 % width=max(mesh.add_x_points);
 % depth=mesh.max_y+mesh.min_y;
 % depth = depth - min(mesh.param_y2);


%
% LOAD JAVA LIBRARY
%
import interp.*;
import edu.mines.jtk.dsp.*;


 %% Constructs a filter with an isotropic Gaussian window of half-width plof
 lof = LocalOrientFilter(plof);

 %% Set half-width of Gaussian derivative filter used to compute gradients
 % lof.setGradientSmoothing(2);

 %% Estimate 2-D structure tensors of training image
 et = lof.applyForTensors(image);

 %% Invert the structure tensors T, ie compute the eigenvalues of the metric (diffusion) tensor field D
 % Before inversion, we have 0 < lv <= lu (if lv=0, then we set lv=eps).
 % After inversion, all eigenvalues are in the range (0,1]. Specifically, after inversion, 0 < au <= av <= 1.
 % If lm is the minimum of the initial eigenvalues lv, 
 % then the parameter p0 is used to compute a0 = pow(lm/lv,p0) = (lm/lv)^p0 
 %  and the parameter p1 is used to compute a1 = pow(lv/lu,p1) = (lv/lu)^p1.
 % Inverted eigenvalues are then
 %  au = a0*a1 = (lm/lv)^p0 * (lv/lu)^p1     and     av = a0 = (lm/lv)^p0,
 % which is consistent with eq.(19) in Zhou et al (2014).
 %
 % p0 emphasizes overall amplitude and p1 emphasizes linearity (ie anisotropy).
 % For amplitude-independent tensors with all eigenvalues av equal to one, set p0 = 0.0.
 % To enhance linearity, set p1 > 1.0. To simply invert (and normalize) these tensors, set p0 = p1 = 1.0.
 et.invertStructure(p0,p1); 
 % JZ: These are the p0, p1 in my paper.
 %     The larger the second parameter, the more anisotropic the tensors,
 %     can be set to 3 or 4, but for geological cross-section image please leave as 1.0.

 %% 2D smoothing directions correspond to eigenvectors of tensors.
 % The direction U corresponds to the largest eigenvalue (perpendicular to linear features),
 % The direction V corresponds to the smallest eigenvalue (parallel to linear features: it is the one of interest for us).
 U = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','U');
 V = javaMethod('valueOf','edu.mines.jtk.dsp.LocalSemblanceFilter$Direction2','V');

 %% Define local smoothing filters for local semblance computation
 %lsf1 = LocalSemblanceFilter(plsf1_1,plsf1_2);   %not used
 lsf2 = LocalSemblanceFilter(plsf2_1,plsf2_2);

 %% Computes local semblance images using local smoothing filters.
 % Local semblance (Hale, 2009) is defined to be a squared smoothed-image divided by a smoothed squared-image,
 % where smoothing is performed by local smoothing filters along the eigenvectors of a structure tensor field.
 % Ref: Hale, D., 2009, Structure-oriented smoothing and semblance, CWP-635
 %s1 = lsf1.semblance(V,et,image);   %not used (would use coherence rather than semblance)
 s2 = lsf2.semblance(U,et,image);

 %% Get eigenvalues of structure tensors
 %FL2JZ: why do you use your own function and not Dave Hale's one?
 %       (my tests suggest it doesn't change anything...)
 %Answer: just because of bugs...
 ev = interp.Sm.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));

 %% Initialize output matrix of structure tensors
 %  (its size is the same as the inverse, coarse, grid,
 %   and btw NOT the same as the training image)
 mesh.et2 = EigenTensors2(m1,m2);

 %% init. scaling matrix (further filled using semblance)
 mesh.s = zeros(m2,m1);


 %% REMIND CONVENTIONS:
 % X = vector of size num_param with all INDICES of grid cells in the x-direction
 % Y =   "    "   "      "       "    "     "    "   "     "   "   "  y-direction
 % mesh.param_x = vector of size num_param with all LOCATIONS of grid cells in the x-direction
 % mesh.param_y =   "    "   "      "       "    "      "     "   "     "   "   "  y-direction
 % XXX = vector of size M1 with LOCATIONS of cells in the x-direction
 % YYY =   "    "   "   M2  "       "     "    "   "   "  y-direction

 %X=ceil((mesh.param_x-min(mesh.param_x))*(m1/width));   %not used anymore
 %Y=ceil((mesh.param_y-min(mesh.param_y))*(m2/depth));
 %%% with topography data, use mesh.param_y2
 %% Y=ceil((mesh.param_y2-min(mesh.param_y2))*(n2/depth));
 XXX=unique(mesh.param_x);
 YYY=unique(mesh.param_y);
 size_XXX=size(XXX)
 size_YYY=size(YYY)


 % LOOP OVER ROWS IN MODEL (y-indices)
 for i2=1:m2

     % LOOP OVER COLUMNS IN MODEL (x-indices)
     for i1=1:m1

         % compute indices (k1,k2) in TI that corresponds to inversion cell (i1,i2)
         k1=ceil(XXX(i1)/input.hstep)+1;
         k2=ceil(YYY(i2)/input.hstep)+1;

         %get eigenvalues at point (k1,k2) in TI
         %(using Jieyi's routine)
         evi = ev(:,k2,k1);

         %this would be using Dave Hale's routine
         %FL2JZ: I could notice in ellipse.m that it gives the same results...
         %evi = et.getEigenvalues(k1-1,k2-1);

         %get the eigenvector U (perpendicular to structures) at point (k1,k2) in TI
         %(be careful to Java indices, starting at 0)
         %FL: we are actually more interested in the eigenvector V (parallel to structures),
         %    but defining the direction of U also define the direction of V
         %    (since U and V are orthogonal)
         %    So mesh.et2 will well contain both U and V.
         %FL2JZ: could you confirm that? YES
         eti = et.getEigenvectorU(k1-1,k2-1);

         %%FL: additional manipulations commented by JZ... (not used)
         %eti = et.getTensor(k1,k2);
         %mesh.et2.setTensor(i1,i2,eti);
         %mesh.et2.invertStructure(1,1);
         %evi = et.getEigenvalues(k1,k2);

         %if abs(a(1)-a(2))<0.99
         %   a(1)=0.25; a(2)=1;
         %   aa(1)=0;   aa(2)=1;
         %end

         %FL2JZ: before you set the 2nd eigenvalue to 1...
         %       I think we should not do that,
         %       because it leads to anisotropic smoothing even in isotropic regions (see attached Fig.).
         %Answer: experimental... not to be used.
         %evi(2)=1;

         % store eigenvalues, eigenvectors and semblance in array mesh.et2 of size [m1,m2] of inversion model
         mesh.et2.setEigenvectorU(i1-1,i2-1,eti);
         mesh.et2.setEigenvalues(i1-1,i2-1,evi);
         mesh.s(i2,i1)=1/(1.0001-s2(k2,k1));
         %mesh.s(i2,i1)=s1(k2,k1);   %FL: alternate definition from JZ... (use coherence rather than semblance)

     end   %i1

 end   %i2

 %% plot ellipses representing structure tensors
 ellipse

end   %end function
