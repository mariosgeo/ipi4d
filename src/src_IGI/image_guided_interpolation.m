%
% This function performs the image-guided interpolation of the recon-
% structed models, using the structure tensors defined in initcm2.m.
% 
% Francois Lavoue and Jieyi Zhou, February 5, 2016

function model_out=image_guided_interpolation(input,mesh,model_in)

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
 % (using the same parameters as in initcm2.m)
 lof = LocalOrientFilter(input.IGI.p_lof1,input.IGI.p_lof2);
 et = lof.applyForTensors(image);
 et.invertStructure(input.IGI.p0,input.IGI.p1);

 % define image size and sampling
 [n2,n1]=size(image);
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
 bg = BlendedGridder2(et,real(model_in),X,Y);
 tmp_real = bg.grid(s1,s2);
 bg = BlendedGridder2(et,imag(model_in),X,Y);
 tmp_imag = bg.grid(s1,s2);
 model_out=tmp_real+1i*tmp_imag;

 % plot interpolated model
 %FL: display has to be checked and improved...
 if input.plot_interpolated_model==1
    figure();
    imagesc(vx,vz,reshape(model_out,[n2,n1]));
    title('INTERPOLATED MODEL');
    colorbar;
 end

