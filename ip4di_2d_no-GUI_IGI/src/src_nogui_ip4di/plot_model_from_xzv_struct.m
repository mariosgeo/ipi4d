%
% Function PLOT_MODEL_FROM_XZV_STRUCT: plot model using the 'mesh' structure.
%
% Derived from function auto_contour.m by M. Karaoulis.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function plot_model_from_xzv_struct(input,mesh,vMX,vMZ,model)

%figure;

% coordinates
% vMX = vector with all x-locations of model parameters
% vMZ =   "     "   all z-loc.      "    "       "
vx=unique(vMX);   % vector of x-locations
vz=unique(vMZ);   %   "    "  z-loc.
nx=length(vx)
nz=length(vz)

% add extra columns on the edges (as mesh does)
%vx=[vx(1)-input.electrode_spacing ; vx ; vx(nx)+input.electrode_spacing];
%model=[model(:,1) model model(:,nx)];


%-- correct for topography
if input.plot_options.plot_topo==1

   % interpolate topo on input x-vector
   ztopo=interp1(mesh.topo_nodes(:,1),mesh.topo_nodes(:,2),vx,'linear');

   for ix=1:nx
       ind=find(vMX==vx(ix));
       if length(ind)==nz
          vMZ(ind)=vMZ(ind)-ztopo(ix);   %FL: WHY DOES "-" WORK AND NOT "+"???
       else
          error(sprintf('Found only %i points at x = %f (ix = %i), there should be nz = %f.',length(ind),vx(ix),ix,nz))
       end
   end

   % new vector of z-locations (does not match the size of model anymore)
   vz=unique(vMZ);
   nz2=length(vz)

%   vix=floor(vMX)+1;
%DataMatrix = accumarray(YourData(:,1:2),YourData(:,3),[],[],FillValue);
end


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


if length(model)>0 && input.plot_options.plot_TI==0
%- plot model
   sz=[];   %max(mesh.dx,mesh.dy);
   hplot=scatter(vMX,-vMZ,sz,model,'filled','s')   %FL: WHY DOES "-" WORK AND NOT "+"???
   tune_figure; %(input,hplot);
   hold on;

elseif length(model)>0 && input.plot_options.plot_TI==1
   sz=[];
   hTI=scatter(vMX,-vMZ,sz,model,'filled','s')
   tune_figure; %(input,hTI);

else
%- plot (x,y)
   sz=[];   %max(mesh.dx,mesh.dy);
   hplot=plot(vMX,-vMZ,'r.');
   tune_figure; %(input,hplot);
   hold on;
end

%%- plot TI with handle "hsc" to make scatterers transparent
%if input.plot_options.plot_TI==1
%   hplot=scatter(vMX,-vMZ,sz,model,'filled','s');
%end

%% tune figure
%tune_figure(input,hplot);

% using imagesc (obsolete)
%model=reshape(model,[nx,nz]);
%model=transpose(model);   %transpose to have it in [nz,nx] format (MK's convention is opposite as SU one)
%imagesc(vx,vz,model);

%elseif input.plot_options.interp==1
%% plot interpolated model
%
%   %ZI=griddata(param_x-neg_x,-param_y,amp,XI,YI,'nearest');
%   ZI=TriScatteredInterp(MX(:),MZ(:),model);
%   ZII = ZI(MX,MZ);
%   contourf(MX,MZ,ZII,17,'EdgeColor','none');
%
%elseif input.plot_options.interp==2
%% plot model after image-guided interpolation
%   ZI=model(:);
%   ZII = ZI(MX,MZ);
%   contourf(MX,MZ,ZII,17,'EdgeColor','none');
%end

end   %end function

