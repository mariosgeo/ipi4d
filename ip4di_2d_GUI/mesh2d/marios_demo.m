clear all

      node = [0.0, 0.0; 1.0,0.0; 1.0,1.0; 0.0,1.0; 2.0,0.0; 2.0,1.0; 3.0,0.0; 3.0,1.0];
      edge = [1,2; 2,3; 3,4; 4,1; 2,5; 5,6; 6,3; 5,7; 7,8; 8,6];
      faces{1} = [1,2,3,4];
      faces{2} = [5,6,7,2];
      faces{3} = [8,9,10,6];
hdata.edgeh=[1,0.05 ;7,1];
 [p,t]=mesh2d(node,edge,hdata);     
%[p,t,fnum]=meshfaces(node,edge,faces,hdata);





p(:,2)=-p(:,2);
node(:,2)=-node(:,2);

%figure('Name','Mesh')
  % plot(x_node_coord,y_node_coord,'b.','markersize',1)
   hold on;
   % Colour mesh for each face
   col = ['b','r','g','k','m'];
   for k = 1:length(faces)
      colk = mod(k,length(col));
      if (colk==0)
         colk = length(col);
      end
      patch('faces',t,'vertices',p,'facecolor','w','edgecolor',col(colk));
   end
  % patch('faces',edge,'vertices',node,'facecolor','none','edgecolor','k')
   % Highlisght low q triangles in debug mode
 %  if options.debug
 %     pc = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3.0;
 %     plot(pc(q<0.5,1),pc(q<0.5,2),'r.')
 %  end
   axis equal;% off;