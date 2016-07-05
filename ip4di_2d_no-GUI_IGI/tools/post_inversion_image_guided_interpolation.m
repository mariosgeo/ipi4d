%
% This program performs the image-guided interpolation
% of the inversion results, using the structure tensors
% extracted for the inversion and stored in structure
% file struct_mesh.mat.
% 
% Francois Lavoue and Jieyi Zhou, February 5, 2016

clear all
close all

%=====   USER PARAMETERS   =====%

% choose to use Matlab or Octave
% (does not work with Octave yet)
matlab_flag=1;   % 1->Matlab, else->Octave

% output file name
file_out_int='interpolated_model.dat';

% add path for plot_model.m
addpath('/home/lavouef/02_SIP/06_SIP_INVERSION/IP4DI_PKG/IP4DI_2D_IGI_no-GUI_v0.1/src/src_nogui_ip4di/');

% path to Dave Hale's JTK and IDH libraries
% (to be changed according to user's system)
LIB_PATH='/home/lavouef/git/';
javaaddpath([LIB_PATH 'idh/bench/build/classes'])
javaaddpath([LIB_PATH 'jtk/build/libs/edu_mines_jtk.jar'])
javaaddpath([LIB_PATH 'jtk/libs/arpack-java.jar'])
javaaddpath([LIB_PATH 'jtk/libs/netlib-java.jar'])
javaaddpath([LIB_PATH 'jtk/libs/gluegen-rt.jar'])
javaaddpath([LIB_PATH 'jtk/libs/jogl-all.jar'])
javaaddpath([LIB_PATH 'jtk/libs/junit.jar'])
javaaddpath([LIB_PATH 'jtk/libs/jythonlib.jar'])

% load structures from run dir.
if matlab_flag==1
   input=importdata('struct_input.mat');
   mesh =importdata('struct_mesh.mat');
else   % Octave syntax
   load('struct_input.mat');
   load('struct_mesh.mat');
   %load('struct_final_lagrn0.01.mat');
end

% option to plot models
input.plot_models=1;
input.plot_options.cmplx_flag=3;   % 3->plot amplitude, 4->plot phase

%=====  END USER PARAMETERS   =====%


% LOAD JAVA LIBRARY
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;

% define size of training image
% /!\ opposite conventions as SU ones
n2=input.n1_TI;   %nz
n1=input.n2_TI;   %nx
vx=[0:n1-1]*input.hstep;
vz=[0:n2-1]*input.hstep;

% read image
fid=fopen(input.training_image,'r');
   image=fread(fid,[n2,n1],'float');
fclose(fid);

% define local smoothing filter
% (using the same parameters input.* as for inversion)
lof = LocalOrientFilter(input.IGI.p_lof1,input.IGI.p_lof2);
et = lof.applyForTensors(image);
et.invertStructure(input.IGI.p0,input.IGI.p1);

% define image sampling
s1 = Sampling(n1,1,0);
s2 = Sampling(n2,1,0);
depth=max(mesh.param_y);
width=max(mesh.param_x);
X=ceil((mesh.param_x-min(mesh.param_x))*(n1/width));
Y=ceil((mesh.param_y-min(mesh.param_y))*(n2/depth));
%(NB: with topography data, use mesh.param_y2)

X(find(X==0))=1;Y(find(Y==0))=1;
X(find(X==n1))=n1-1;Y(find(Y==n2))=n2-1;
X(find(X==n1+1))=n1-1;Y(find(Y==n2+1))=n2-1;

% perform structure-guided interpolation
bg = BlendedGridder2(et,real(mesh.res_param1),X,Y);
interpolated_model_real = bg.grid(s1,s2);   %real part
bg = BlendedGridder2(et,imag(mesh.res_param1),X,Y);
interpolated_model_imag = bg.grid(s1,s2);   %imag. part
interpolated_model = interpolated_model_real + 1i*interpolated_model_imag;
interpolated_model = reshape(interpolated_model,[n2,n1]);

% save interpolated model
fid=fopen(file_out_int,'w');
for i1=1:n1
    for i2=1:n2
        if input.cmplx_format==1
        %save real and imaginary part of interpolated resistivity
           fprintf( fid,'%f %f %f %f\n',...
             vx(i1),vz(i2),real(interpolated_model(i2,i1)),imag(interpolated_model(i2,i1)) );
        elseif input.cmplx_format==2
        %save amplitude and phase of interpolated resistivity
           fprintf( fid,'%f %f %f %f\n',...
             vx(i1),vz(i2),abs(interpolated_model(i2,i1)),1000*atan2(imag(interpolated_model(i2,i1)),real(interpolated_model(i2,i1))) );
        end
    end   %i2
end   %i1
fclose(fid);

% plot models
if input.plot_models==1
   % plot inversion result
   figure()
   input.plot_options.label_title='INVERSION RESULT';
   plot_model_on_forward_mesh(input,mesh,mesh.res_param1);

   % plot interpolated model
   figure()
   imagesc(vx,vz,interpolated_model);
   colorbar
end

