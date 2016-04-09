%
% Function INITCM2: build prior model covariance matrix with
%                   structured-based local directional Laplacian filters
%
% Author: Jieyi Zhou, Colorado School of Mines, 2014
%
% Modified by F. Lavoue', Colorado School of Mines, September 2015
% in order to:
% 1) consider TI and grid of different sizes.
% 2) comment the use of Dave Hale's and Jieyi Zhou's Java routines
%    (comments are based on Dave Hale's documentation at http://dhale.github.io/jtk/api)
% 3) build the covariance matrix explicitly,
%    which enables a) to compute the resolution matrix with the contribution of the prior information,
%                  b) to apply spatially-varying Lagrangian multipliers via ACB (Yi et al, 2003).

function [mesh]=initcm2(input,mesh,cx,cz)

disp(' ')
disp('------------------------------------------------------')
disp('                   ENTER INITCM2                      ')
disp(' Define structure-based directional Laplacian filters ')
disp('------------------------------------------------------')
disp(' ')

 %%=== TUNING PARAMETERS ===%%

 %factor = 2;   %JZ: this is just for how many image samples you will use to calculate a tensor,
 %                   the larger the factors, the less samples are used
 %FL: not used anymore

 %mesh.scale=1;   %FL: should probably not be used anymore (ie let to 1),
                  %    because it is redundant with the spatially varying Lagrangian distribution defined by ACB.
                  % (now defined in 'inversion_parameters.m')

 p_lof1=input.IGI.p_lof1;   % half-widths of isotropic Gaussian window for LocalOrientFilter
 p_lof2=input.IGI.p_lof2;

 %p_grad=input.IGI.p_grad;   % half-width of Gaussian derivative filter used to compute gradients in LocalOrientFilter
                             % for the calculation of the structure tensor T
 %FL: should not be modified (according to Elias Arias' experience...)

 p0=input.IGI.p0;    %JZ: These are the p0, p1 in my paper.
 p1=input.IGI.p1;    %    The larger the second parameter, the more anisotropic the tensors,
                     %    can be set to 3 or 4, but for geological cross-section image please leave as 1.0.
 %DH's doc for invertStructure:
 % p0 emphasizes overall amplitude and p1 emphasizes linearity (ie anisotropy).
 % For amplitude-independent tensors with all eigenvalues av equal to one, set p0 = 0.0.
 % To enhance linearity, set p1 > 1.0. To simply invert (and normalize) these tensors, set p0 = p1 = 1.0.

 stencil='D33';   % stencil used for computing derivatives in LocalDiffusionKernel (D22, D33, D24, D71 or D91)
 %NB: hard-coded in Java routine 'java-src/IGICovariance/pkgApplyGDG.java' for the moment, but should be passed as argument on the long term.

 %% Parameters of local smoothing filters (lsf) for local semblance computation
 p_lsf1_1=input.IGI.p_lsf1_1;   % lsf1 is related to coherence (not used)
 p_lsf1_2=input.IGI.p_lsf1_2;

 p_lsf2_1=input.IGI.p_lsf2_1;   % half-width of 1st smoothing filter for lsf2 (related to semblance, used)
 p_lsf2_2=input.IGI.p_lsf2_2;   % half-width of 2nd smoothing filter

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
 mesh.m1=m1;   %nx
 mesh.m2=m2;   %nz

 if input.debug==1
 % compare domain sizes of TI and inversion model
 % (a slight difference is normal due to the location
 %  of element centers reported in param_x,param_y)
    depth=max(vy_TI)
    width=max(vx_TI)
    depth=max(mesh.param_y)
    width=max(mesh.param_x)
 end


%
% LOAD JAVA LIBRARY
%
import IGICovariance.*;
import edu.mines.jtk.dsp.*;


 %% Constructs a filter with an isotropic Gaussian window of half-width p_lof
 %% or with an anisotropic Gaussian window of half-widths p_lof1,p_lof2
 lof = LocalOrientFilter(p_lof1,p_lof2);
 %%NB: the half-widths p_lofs relate to a second smoothing applied after building
 %%    the tensor T (K_rho in Weckert, )
 %%FL: according to Elias Arias' experience on dip filtering, an anisotropic Gaussian
 %%    window is suitable if the training image contains anisotropic structures.
 %%    If the structures display more variations in direction 1 (vertical), we may then
 %%    use a larger p_lof1 than p_lof2 (e.g. [p_lof1,p_lof1/2] or [p_lof1,p_lof1/3]).

 %% Set half-width of Gaussian derivative filter used to compute gradients
 %% This is a first smoothing applied when computing the derivatives using Gaussian windows (K_sigma in Weckert)
 %%FL: should not be changed according to Elias Arias' experience.
 %lof.setGradientSmoothing(p_grad);

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
 lsf1 = LocalSemblanceFilter(p_lsf1_1,p_lsf1_2);   %not used
 lsf2 = LocalSemblanceFilter(p_lsf2_1,p_lsf2_2);

 %% Computes local semblance images using local smoothing filters.
 % Local semblance (Hale, 2009) is defined to be a squared smoothed-image divided by a smoothed squared-image,
 % where smoothing is performed by local smoothing filters along the eigenvectors of a structure tensor field.
 % Ref: Hale, D., 2009, Structure-oriented smoothing and semblance, CWP-635
 s1 = lsf1.semblance(V,et,image);   %coherence (not used)
 s2 = lsf2.semblance(U,et,image);   %semblance (used to make weak structures appear quasi-isotropic)

 %% Get eigenvalues of structure tensors
 %FL: this uses JZ's routine instead of DH's one to get eigenvalues at all locations
 %    (instead of local eigenvalues only)
 ev = IGICovariance.pkgApplyGDG.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));

 %% Initialize Java array of structure tensors
 %  (its size is the same as the inverse, coarse, grid,
 %   and NOT the same as the training image)
 % (be careful to Java indices, in reverse order compared to Matlab ones)
 mesh.et2 = EigenTensors2(m1,m2);

 %% init. scaling matrix (further filled using semblance)
 mesh.s = zeros(m2,m1);

 %% REMIND CONVENTIONS:
 % mesh.param_x = vector of size num_param with all LOCATIONS of grid cells in the x-direction
 % mesh.param_y =   "    "   "      "       "    "      "     "   "     "   "   "  y-direction
 % XX = vector of size M1 with LOCATIONS of cells in the x-direction
 % YY =   "    "   "   M2  "       "     "    "   "   "  y-direction
 XX=unique(mesh.param_x);
 YY=unique(mesh.param_y);


 % LOOP OVER ROWS IN MODEL (z-indices)
 for i2=1:m2

     % LOOP OVER COLUMNS IN MODEL (x-indices)
     for i1=1:m1

         % find indices (k1,k2) in TI that corresponds to inversion cell (i1,i2)
         k1=ceil(XX(i1)/input.hstep)+1;
         k2=ceil(YY(i2)/input.hstep)+1;

         % get eigenvalues at point (k1,k2) in TI
         % use JZ's routine
         evi = ev(:,k2,k1);
         %evi = et.getEigenvalues(k1-1,k2-1);   %FL: DH's routine does not work anymore...?

         %FL TO DO
         % scale eigenvalues using semblance
         % such as to make weak structures quasi-isotropic
         %evi(1) = evi(1)*??;
         %evi(2) = evi(2)*??; 

         % get eigenvectors at point (k1,k2) in TI
         % (be careful to Java indices, starting at 0,
         %  and in reverse order compared to Matlab ones)
         ui = et.getEigenvectorU(k1-1,k2-1);
         vi = et.getEigenvectorV(k1-1,k2-1);

         % keep horizontal and vertical components of U and V
         % to build Cm later on
         Ux(i2,i1)=evi(1)*ui(1);
         Uz(i2,i1)=evi(1)*ui(2);
         Vx(i2,i1)=evi(2)*vi(1);
         Vz(i2,i1)=evi(2)*vi(2);

         % store eigenvalues, eigenvectors and semblance in array mesh.et2 of size [m1,m2] of inversion model
         %NB: storing the eigenvector U is sufficient since V is orthogonal
         mesh.et2.setEigenvectorU(i1-1,i2-1,ui);
         mesh.et2.setEigenvalues(i1-1,i2-1,evi);
         mesh.s(i2,i1)=1/(1.0001-s2(k2,k1));     %FL: not sure it is the right formulation, seems quite unstable...
         %mesh.s(i2,i1)=s1(k2,k1);   %FL: alternate definition from JZ... (use coherence rather than semblance)
     end   %i1
 end   %i2


 %% BUILD COVARIANCE MATRIX
 %  as Cm^-1 =        Wu'Wu        +        Wv'Wv
 %           = Ux*Wx'Wx + Uz*Wz'Wz + Vx*Wx'Wx + Vz*Wz'Wz
 mesh.ctc = zeros(mesh.num_param,mesh.num_param);   % init.
 Ux=Ux(:); Uz=Uz(:); Vx=Vx(:); Vz=Vz(:);   % recast Ux,Vx as vectors of size num_param

 % directional derivatives
 WUx = diag(Ux)*cx;
 WUz = diag(Uz)*cz;
 WuWu= WUx'*WUx + WUz'*WUz;

 WVx = diag(Vx)*cx;
 WVz = diag(Vz)*cz;
 WvWv= WVx'*WVx + WVz'*WVz;

 % covariance matrix
 mesh.ctc = WuWu+WvWv;

 %% plot ellipses representing structure tensors
 ellipse;

 %% OUTPUT TMP
 mesh.s1=s1;
 mesh.s2=s2;

end   %end function
