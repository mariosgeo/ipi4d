%
% PROGRAM TO COMPUTE LOCAL DIP FILTERS FROM AN IMAGE
%

clear all
close all

%% PATH BELOW IS ALREADY IN STATIC PATH
%javaaddpath('/Users/lavoue/git/idh/bench/build/classes')
%javaaddpath('/Users/lavoue/git/jtk/build/libs/edu_mines_jtk.jar')
%javaaddpath('/Users/lavoue/git/jtk/libs/arpack-java.jar')
%javaaddpath('/Users/lavoue/git/jtk/libs/netlib-java.jar')
%javaaddpath('/Users/lavoue/git/jtk/libs/gluegen-rt.jar')
%javaaddpath('/Users/lavoue/git/jtk/libs/jogl-all.jar')
%javaaddpath('/Users/lavoue/git/jtk/libs/junit.jar')
%javaaddpath('/Users/lavoue/git/jtk/libs/jythonlib.jar')


%% TRAINING IMAGE
%input.training_image='sandbox_permittivity_v1_125x300_h0.004.bin'
input.training_image='sandbox_permittivity_v1_51x121_h0.01.bin'
input.training_image='sandbox_permittivity_v2_no-box_air_125x300_h0.004.bin'
input.training_image='rtm_sandbox_v2_no-box_air_125x300_h0.004.bin'
input.training_image='gradient_init_sandbox_v2_125x300_h0.004.bin'

input.n1_TI=125
input.n2_TI=300
input.hstep=0.004

%air layer (not used)
%h_air=0.1

%under-sampling rate to avoid memory issues
mesh.usamp=1;

%J. Zhou: There are two ways you can do the image-guided inversion, first is the same as my paper, use the directions of the eigenvectors and do four-direction smoothing, by calling the function initcm1 before inversion to initialize the smoothing matrix. The input to this function includes cx, cy, cd1, cd2 these four first-order derivative matrices with -1 and 1 pairs (see test_jtk_initcm1.m).


%%
%% USE INITCM1.m
%%
% define mesh for using initCm1
% J. Zhou: mesh is just a structure carrying parameters around
% please pass your own parameter and replace them
% mesh.m2 = number of rows of model matrix
% mesh.m1 = number of columns of model matrix
% depth = total depth of model section
% width = total width of model section
% mesh.param_x is the x-axis location vector of each cell of model, row after row
% mesh.param_y is the y-axis location (depth) vector of each cell of model, row after row
% mesh.param_y2 is the depth vector of each cell of model with topography, row after row
% 
% There are some parameters that you can tune, please see the comments.
% 

% /!\ opposite conventions...
mesh.m2=n1
mesh.m1=n2

depth=(n1-1)*hstep
width=(n2-1)*hstep

mesh.param_x=[];
mesh.param_y=[];
%mesh.param_y2=[];   %used instead of param_y in case of topo

for i2=1:mesh.usamp:n2
    tmp_x1=(i2-1)*hstep*ones(n1,1);
    tmp_y1=hstep*[0:n1-1]';

    %under-sample for memory issues
    tmp_x2=tmp_x1(1:mesh.usamp:end);
    tmp_y2=tmp_y1(1:mesh.usamp:end);

    mesh.param_x=[mesh.param_x ; tmp_x2];
    mesh.param_y=[mesh.param_y ; tmp_y2];
    %mesh.param_y2=[mesh.param_x ; hstep*[0:n1-1]'-h_air];
end
clear tmp_x1 tmp_y1 tmp_x2 tmp_y2


%% use initcm1 function 
%[mesh,ctc,s1,s2]=initcm1(mesh,cx,cy,cd1,cd2);


%%only if needed because large
%file_out='tmp_ctc.bin'
%fid=fopen(file_out,'w');
%  fwrite(fid,full(ctc),'real*4');
%fclose(fid);


% TEST2: solve trivial inverse problem min Cm^-1 (m-m_prior)


%%
%% USE INITCM2.m
%%
% The second way is to use a local Laplacian operator, so 
%   Cm = C' D C, 
% where D is the tensor field, C is the gradient operator. This can be done be calling the function initcm2 before inversion to initialize the local Laplacian operator. Then during the inversion, instead of using 
%   dX = (J'J + beta*Cm)^(-1) * J'*Cd*[G(X)-d)], 
% do this: 
%   dX = reshape(interp.Sm.getLaplacianI(mesh.et2,reshape(J'*Cd*[G(X) - d)],mesh.m1,mesh.m2)',(J'J)',scale,mesh.s)',mesh.num_param,1);
% Here scale is how you want to scale the tensors, the larger the more emphasize on structure mesh.s is using structur-oriented semblance s1 to scale the tensors, if don't wanna use can just put [] instead of mesh.s
% There are some parameters that you can tune and they are in the comments.


% use initcm1 function 
[mesh,s2]=initcm2(mesh);
whos


%% OUTPUT semblance and scaling matrices
%file_out='tmp_s1.bin'
%fid=fopen(file_out,'w');
%  fwrite(fid,full(s1),'real*4');
%fclose(fid);

file_out='tmp_s2.bin'
fid=fopen(file_out,'w');
  fwrite(fid,full(s2),'real*4');
fclose(fid);

file_out='tmp_meshs.bin'
fid=fopen(file_out,'w');
  fwrite(fid,single(full(mesh.s)),'real*4');
fclose(fid);

