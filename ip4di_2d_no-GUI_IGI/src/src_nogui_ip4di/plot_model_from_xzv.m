%
% Function PLOT_MODEL: simply plot model when mesh is not known
%
% Derived from function auto_contour.m by M. Karaoulis.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function plot_model_from_xzv(input,vMX,vMZ,model)

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


if input.plot_options.cmplx_flag==1
% plot real part
   model=real(model);

elseif input.plot_options.cmplx_flag==2
% plot imag. part
   model=imag(model);

elseif input.plot_options.cmplx_flag==3
% plot amplitude
   model=abs(model);

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
   disp('POST-INVERSION IMAGE-GUIDED INTERPOLATION NOT IMPLEMENTED YET');
   return;
end


% tune figure
colormap(input.plot_options.cmap);
hleg=colorbar();
title(input.plot_options.label_title);
axis ('equal');
xlabel('x (m)');
ylabel('z (m)');
xlim([input.plot_options.x_min input.plot_options.x_max]);
ylim([input.plot_options.z_min input.plot_options.z_max]);

% tune figure according to plotted component
if input.plot_options.cmplx_flag==1
   ylabel(hleg,'Resistivity, real part (\Omega.m)')

elseif input.plot_options.cmplx_flag==2
   ylabel(hleg,'Resistivity, imag. part (\Omega.m)')

elseif input.plot_options.cmplx_flag==3
   ylabel(hleg,'Resistivity, amplitude (\Omega.m)')
   input.plot_options.caxis=input.plot_options.caxis_amp;
   input.plot_options.axis_tics=input.plot_options.axis_tics_amp;

elseif input.plot_options.cmplx_flag==4
   ylabel(hleg,'Resistivity, phase (mrad)')
   input.plot_options.caxis=input.plot_options.caxis_phi;
   input.plot_options.axis_tics=input.plot_options.axis_tics_phi;
end

if(numel(input.plot_options.caxis)>0); caxis(input.plot_options.caxis); end
if(numel(input.plot_options.axis_tics)>0); set(hleg,'YTick',input.plot_options.axis_tics); end

end   %end function

