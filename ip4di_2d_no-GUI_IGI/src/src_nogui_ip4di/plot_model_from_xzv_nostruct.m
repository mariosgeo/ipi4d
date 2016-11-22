%
% Function PLOT_MODEL_FROM_XZV_NOSTRUCT: plot model when the 'mesh'
%                                        structure is not known.
%
% Derived from function auto_contour.m by M. Karaoulis.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function plot_model_from_xzv_nostruct(input,vMX,vMZ,model)

figure;

%-----------------------Plots----------------------------------------------

% coordinates
% vMX = vector with all x-locations of model parameters
% vMZ =   "     "   all z-loc.      "    "       "
vx=unique(vMX);   % vector of x-locations
vz=unique(vMZ);   %   "    "  z-loc.
nx=length(vx)
nz=length(vz)

model=reshape(model,[nx,nz]);
model=transpose(model);   %transpose to have it in [nz,nx] format (MK's convention is opposite as usual one)

% add extra columns on the edges (as mesh does)
vx=[vx(1)-input.electrode_spacing ; vx ; vx(nx)+input.electrode_spacing];
model=[model(:,1) model model(:,nx)];

% create mesh
[MX,MZ]=meshgrid(vx,vz);


if input.plot_options.cmplx_flag==0 || input.plot_options.cmplx_flag==1
% plot real part
   if input.plot_options.plot_log==0
      model=real(model);
   else
      model=real(log10(model));
   end

elseif input.plot_options.cmplx_flag==2
% plot imag. part
   if input.plot_options.plot_log==0
      model=imag(model);
   else
      model=imag(log10(model));
   end

elseif input.plot_options.cmplx_flag==3
% plot amplitude
   if input.plot_options.plot_log==0
      model=abs(model);
   else
      model=abs(log10(model));
   end

elseif input.plot_options.cmplx_flag==4
% plot phase
   model=1000*atan(imag(model)./real(model));
end


if input.plot_options.interp==0
% plot model with inversion discretization
   imagesc(vx,vz,model);

elseif input.plot_options.interp==1
% plot interpolated model

   %ZI=griddata(param_x-neg_x,-param_y,amp,XI,YI,'nearest');
   ZI=TriScatteredInterp(MX(:),MZ(:),model(:));
   ZII = ZI(MX,MZ);
   contourf(MX,MZ,ZII,17,'EdgeColor','none');

elseif input.plot_options.interp==2
% plot model after image-guided inerpolation
   disp('POST-INVERSION IMAGE-GUIDED INTERPOLATION NEEDS THE MESH STRUCTURE');
   disp('USE FUNCTION plot_model_from_xzv_struct(input,mesh,vx,vz,model)');
   return;
end


% tune figure
tune_figure;


end   %end function

