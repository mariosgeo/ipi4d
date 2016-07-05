% reverse z-axis to negative values
% (just for plotting...)
%p(:,2)=-p(:,2);
%node(:,2)=-node(:,2);

%figure('Name','Mesh')
if input.plot_mesh==1
   figure
   %plot(mesh.x_node_coord,mesh.y_node_coord,'b.','markersize',1)
   hold on;

   % Colour mesh for each face
   col = ['b','r','g','k','m'];
   for k = 1:length(faces)
      colk = mod(k,length(col));
      if (colk==0)
         colk = length(col);
      end
       patch('faces',t(fnum==k,:),'vertices',p,'facecolor','w','edgecolor',col(colk));
   end
   patch('faces',edge,'vertices',node,'facecolor','none','edgecolor','k')

   % Highlisght low q triangles in debug mode
   %if options.debug
   %   pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3.0;
   %   plot(pc(q<0.5,1),pc(q<0.5,2),'r.')
   %end

   % tune figure
   title('MESH')
   set(gca,'YDir','reverse');

   axis equal;% off;

   xlim([min(mesh.x_node_coord),max(mesh.x_node_coord)])
   ylim([min(mesh.y_node_coord),max(mesh.y_node_coord)])

   %set(gca,'XTick',[min(mesh.x_node_coord):0.1:max(mesh.x_node_coord)]);
   %set(gca,'YTick',[min(mesh.y_node_coord):0.1:max(mesh.y_node_coord)]);

   xlabel('x (m)')
   ylabel('z (m)')

end   %end if plot_mesh
