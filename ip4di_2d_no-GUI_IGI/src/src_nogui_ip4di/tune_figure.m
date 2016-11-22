%
% SUBROUTINE TUNE_FIGURE: a unique routine to tune all figures consistently.
%
% Using the 'input' Matlab structure and the input.plot_options variables.
% NB: it does not need the 'mesh' structure.
%
% Francois Lavoue, 6 Oct. 2016
%

%function tune_figure(input,hfig)
%FL: not defined as a function anymore but as a simple routine, enabling
%    to pass any variables from the calling program transparently.

figure(gcf)
%set(gcf, 'Color', 'None')

%-- tune figure according to plotted component
% default is amplitude
input.plot_options.caxis=input.plot_options.caxis_amp;
input.plot_options.axis_tics=input.plot_options.axis_tics_amp;

if input.plot_options.plot_TI==0

   % add colorbar
   hleg=colorbar;

   if input.plot_options.cmplx_flag==0
   %- DC case: just resistivity
      ylabel(hleg,'Resistivity (\Omega.m)')

   elseif input.plot_options.cmplx_flag==1
   %- real part
      ylabel(hleg,'Resistivity, real part (\Omega.m)')

   elseif input.plot_options.cmplx_flag==2 && input.sip_flag==1
   %- imag part
      ylabel(hleg,'Resistivity, imag. part (\Omega.m)')

   elseif input.plot_options.cmplx_flag==3
   %- amplitude
      ylabel(hleg,'Resistivity, amplitude (\Omega.m)')

   elseif input.plot_options.cmplx_flag==4 && input.sip_flag==1
   %- phase
      ylabel(hleg,'Resistivity, phase (mrad)')
      input.plot_options.caxis=input.plot_options.caxis_phi;
      input.plot_options.axis_tics=input.plot_options.axis_tics_phi;
   end

   if(numel(input.plot_options.axis_tics)>0); set(hleg,'YTick',input.plot_options.axis_tics); end
end

%-- modify color scale
if(numel(input.plot_options.cmap)>0); colormap(input.plot_options.cmap); end
if(numel(input.plot_options.caxis)>0); caxis(input.plot_options.caxis); end

%-- axis
axis equal;
xlabel('x (m)')
ylabel('z (m)')
if(input.plot_options.reverse_yaxis==1); set(gca,'YDir','reverse'); end
if(numel(input.plot_options.label_title)>0); title(input.plot_options.label_title); end

%-- adjust figure size
xlim([input.plot_options.x_min input.plot_options.x_max]);
ylim([input.plot_options.z_min input.plot_options.z_max]);

%-- log scale
if input.plot_options.plot_log==1
   caxis(log10(input.plot_options.caxis))
   set(hleg,'YTick',log10(input.plot_options.axis_tics));
   set(hleg,'YTickLabel',input.plot_options.axis_tics);
end

%-- adjust transparency if plotting TI
if input.plot_options.plot_TI==1

   %% Opt. 2
   %http://undocumentedmatlab.com/blog/plot-markers-transparency-and-color-gradient
   alpha=0.1
   drawnow
   hMarkers = hTI.MarkerHandle;  % a matlab.graphics.primitive.world.Marker object
   hMarkers.get
   %hMarkers.FaceColorData = uint8(255*[1;0;0;0.3]);  % Alpha=0.3 => 70% transparent red
   hMarkers.FaceColorData(4,:)=uint8(255*alpha);
end

%if numel(input.plot_options.x_min)>0 && numel(input.plot_options.x_max)>0
%   xlim([input.plot_options.x_min input.plot_options.x_max]);
%else
%   xlim([min(mesh.tmp_param(:,3)) max(mesh.tmp_param(:,4))]);
%end
%if numel(input.plot_options.z_min)>0 && numel(input.plot_options.z_max)>0
%   ylim([input.plot_options.z_min input.plot_options.z_max]);
%elseif input.topo_flag==0
%   ylim([min(mesh.tmp_param(:,5)) max(mesh.tmp_param(:,6))]);
%end

