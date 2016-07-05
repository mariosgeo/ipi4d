%
% FUNCTION PLOT_MODEL_ON_FORWARD_MESH
%
% Plots a given model on the cartesian grid encompassing
% the finite-element forward mesh
%
% NB: the variables input.plot_options.cmplx_flag and input.plot_options.label_title
%     must be defined before calling the function.
%
% Francois Lavoue', Colorado School of Mines, October 15, 2015

function plot_model_on_forward_mesh(input,mesh,model)

figure
hold on

for i=1:mesh.num_param
    tmp=mesh.faces{i};
    tmp2=mesh.edge(tmp,:);
    tmp3=mesh.node(tmp2,:);
    tmp4=unique(tmp3,'rows');

    a=find(tmp4(:,2)==max(tmp4(:,2))); % up

    b=find(tmp4(:,2)==min(tmp4(:,2))); % up
    b=b(end:-1:1);

    c=find(tmp4(:,1)==min(tmp4(:,1))); % left
    cc=setdiff(c,a);
    c=setdiff(cc,b);

    d=find(tmp4(:,1)==max(tmp4(:,1))); % right
    dd=setdiff(d,a);
    d=setdiff(dd,b);

    faces=[a;d;b;c];

    if input.plot_options.cmplx_flag==1
       patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',real(model(i))','edgecolor','k');
    elseif input.plot_options.cmplx_flag==2
       patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',imag(model(i))','edgecolor','k');
    elseif input.plot_options.cmplx_flag==3
       patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',abs(model(i))','edgecolor','k');
    elseif input.plot_options.cmplx_flag==4
       patch('faces',faces','vertices',tmp4,'facecolor','flat','FaceVertexCData',atan(imag(model(i))/real(model(i)))*1000','edgecolor','k');
    end

end   %end for i=1:mesh.num_param

% plot electrode locations
plot(mesh.orig_probe_x,mesh.orig_probe_z,'o');

axis equal;
hleg=colorbar;

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

if(numel(input.plot_options.cmap)>0); colormap(input.plot_options.cmap); end
if(numel(input.plot_options.caxis)>0); caxis(input.plot_options.caxis); end
if(numel(input.plot_options.axis_tics)>0); set(hleg,'YTick',input.plot_options.axis_tics); end
if(numel(input.plot_options.label_title)>0); title(input.plot_options.label_title); end

xlabel('x (m)')
ylabel('z (m)')

set(gca,'YDir','reverse');

xlim([min(mesh.tmp_param(:,3)) max(mesh.tmp_param(:,4))]);
ylim([min(mesh.tmp_param(:,5)) max(mesh.tmp_param(:,6))]);
%ylim([-max(mesh.tmp_param(:,6)) -min(mesh.tmp_param(:,5))]);

%set(gca,'XTick',[min(mesh.tmp_param(:,3)):0.1:max(mesh.tmp_param(:,4))]);
%set(gca,'YTick',[min(mesh.tmp_param(:,5)):0.1:max(mesh.tmp_param(:,6))]);

