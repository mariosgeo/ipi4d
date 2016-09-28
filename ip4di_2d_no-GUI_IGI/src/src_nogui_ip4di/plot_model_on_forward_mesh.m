%
% FUNCTION PLOT_MODEL_ON_FORWARD_MESH
%
% Plots a given model on the cartesian grid encompassing
% the finite-element forward mesh
%
% NB: the variables input.plot_options.cmplx_flag and input.plot_options.label_title
%     must be defined before calling the function.
%
% Francois Lavoue, Colorado School of Mines, October 15, 2015
%
% Modified FL, August 2016
%  -> take topographic distortion of parameter cells into account.

function plot_model_on_forward_mesh(input,mesh,model)

figure
hold on

%-- color (red, green, black)
col=[1 0 0 ; 0 1 0 ; 1 1 1];


%-- loop over parameter cells
for i=1:mesh.num_param

    tmp=mesh.faces{i};         % list of node/edge indices for parameter cell i
    tmp2=mesh.edge(tmp,:);     % edge indices
    tmp3=mesh.node(tmp2,:);    % nodes coordinates
    tmp4=unique(tmp3,'rows');  %   "       "

    x1=min(tmp4(:,1));         % x cell boundaries
    x2=max(tmp4(:,1));
    l1=find(tmp4(:,1)==x1);    % and their indices within tmp4
    l2=find(tmp4(:,1)==x2);

    % y cell boundaries
    y1min=tmp4(l1(find(tmp4(l1,2)==min(tmp4(l1,2)))),2);
    y1max=tmp4(l1(find(tmp4(l1,2)==max(tmp4(l1,2)))),2);
    y2min=tmp4(l2(find(tmp4(l2,2)==min(tmp4(l2,2)))),2);
    y2max=tmp4(l2(find(tmp4(l2,2)==max(tmp4(l2,2)))),2);

    % cell vertices
    p1=find(tmp4(:,1)==x1 & tmp4(:,2)==y1max);   % top left
    p2=find(tmp4(:,1)==x2 & tmp4(:,2)==y2max);   % top right
    p3=find(tmp4(:,1)==x2 & tmp4(:,2)==y2min);   % bottom right
    p4=find(tmp4(:,1)==x1 & tmp4(:,2)==y1min);   % bottom left

    faces=[p1;p2;p3;p4];
    icol=1+mod(i,3);

    if input.plot_options.plot_log==0
    %- plot resistivity

       if input.plot_options.cmplx_flag==0 || input.plot_options.cmplx_flag==1
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',real(model(i))','edgecolor','k');
       elseif input.plot_options.cmplx_flag==2
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',imag(model(i))','edgecolor','k');
       elseif input.plot_options.cmplx_flag==3
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',abs(model(i))','edgecolor','k');
          %patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',abs(model(i))','edgecolor',col(icol,:));
       elseif input.plot_options.cmplx_flag==4
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',atan(imag(model(i))/real(model(i)))*1000','edgecolor','k');
       end

    elseif input.plot_options.plot_log==1
    %- plot log(resistivity)

       if input.plot_options.cmplx_flag==0 || input.plot_options.cmplx_flag==1
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',real(log10(model(i)))','edgecolor','k');
       elseif input.plot_options.cmplx_flag==2
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',imag(log10(model(i)))','edgecolor','k');
       elseif input.plot_options.cmplx_flag==3
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',abs(log10(model(i)))','edgecolor','k');
          %patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',abs(model(i))','edgecolor',col(icol,:));
       elseif input.plot_options.cmplx_flag==4
          patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',atan(imag(model(i))/real(model(i)))*1000','edgecolor','k');
       end
    end

end   %end for i=1:mesh.num_param


%- plot electrode locations
plot(mesh.orig_probe_x,mesh.orig_probe_z,'o');

if input.debug>0
%- plot node locations
   scatter(mesh.node(:,1),mesh.node(:,2),'.','r');
end

axis equal;
hleg=colorbar;

%-- tune figure according to plotted component
% default is amplitude
input.plot_options.caxis=input.plot_options.caxis_amp;
input.plot_options.axis_tics=input.plot_options.axis_tics_amp;

if input.plot_options.cmplx_flag==0
%- DC case: just resistivity
   ylabel(hleg,'Resistivity (\Omega.m)')

elseif input.plot_options.cmplx_flag==1
   ylabel(hleg,'Resistivity, real part (\Omega.m)')

elseif input.plot_options.cmplx_flag==2 && input.sip_flag==1
   ylabel(hleg,'Resistivity, imag. part (\Omega.m)')

elseif input.plot_options.cmplx_flag==3
   ylabel(hleg,'Resistivity, amplitude (\Omega.m)')

elseif input.plot_options.cmplx_flag==4 && input.sip_flag==1
   ylabel(hleg,'Resistivity, phase (mrad)')
   input.plot_options.caxis=input.plot_options.caxis_phi;
   input.plot_options.axis_tics=input.plot_options.axis_tics_phi;
end

if(numel(input.plot_options.cmap)>0); colormap(input.plot_options.cmap); end
if(numel(input.plot_options.caxis)>0); caxis(input.plot_options.caxis); end
if(numel(input.plot_options.axis_tics)>0); set(hleg,'YTick',input.plot_options.axis_tics); end
if(numel(input.plot_options.label_title)>0); title(input.plot_options.label_title); end

xlabel('x (m)')
ylabel('z (m)')

set(gca,'YDir','reverse');

%
xlim([min(mesh.tmp_param(:,3)) max(mesh.tmp_param(:,4))]);
if numel(input.plot_options.ymin)>0 && numel(input.plot_options.ymax)>0
   ylim([input.plot_options.ymin input.plot_options.ymax]);
elseif input.topo_flag==0
   ylim([min(mesh.tmp_param(:,5)) max(mesh.tmp_param(:,6))]);
end

%ylim([min(mesh.tmp_param(:,5))-min(mesh.anaglyfo_data(:,2)) max(mesh.tmp_param(:,6))-max(mesh.anaglyfo_data(:,2))]);
%ylim([-max(mesh.tmp_param(:,6)) -min(mesh.tmp_param(:,5))]);

%set(gca,'XTick',[min(mesh.tmp_param(:,3)):0.1:max(mesh.tmp_param(:,4))]);
%set(gca,'YTick',[min(mesh.tmp_param(:,5)):0.1:max(mesh.tmp_param(:,6))]);

%-- log scale
if input.plot_options.plot_log==1
   caxis(log10(input.plot_options.caxis))
   set(hleg,'YTick',log10(input.plot_options.axis_tics));
   set(hleg,'YTickLabel',input.plot_options.axis_tics);
end

%FL DEBUG
if input.debug>0
   disp_axis_tics=input.plot_options.axis_tics
   disp_caxis=input.plot_options.caxis
   min_rho=min(model)
   max_rho=max(model)
end
