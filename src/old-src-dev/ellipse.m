r = 1; %parameter to scale the size of ellipse                         
nt = 101; %discretize theta for plotting
dt = 2.0*pi/(nt-1);
ft = 0;
x1=zeros(nt,1);
x2=zeros(nt,1);
c=0;
YYY=unique(mesh.param_y);
% ev = interp.Sm.getEigenvalues(mesh.et2,zeros(m2,m1),zeros(m2,m1));
%     for i=1:mesh.num_param
%       i1 = mod(i,m1);if(i1==0) i1=m1-1; else i1=i1-1; end
%     i2 = ceil(i/m1)-1; 
%         v = mesh.et2.getEigenvectorV(i1,i2);
%         s = mesh.et2.getEigenvalues(i1,i2);
% % s = ev(:,i2,i1);
% 
%         v1 = v(1); v2 = v(2);
%         du = s(1); dv = s(2);
% % if abs(s(1)-s(2))<1
% % %     du=0.25; dv=1;
% %     v1=1;v2=0;
% % end
%     i1=mesh.param_x(i1+1);
%           i2=YYY(i2+1);
%         a = r*sqrt(dv); b = r*sqrt(du);
%         for it=1:nt
%           c=c+1;
%           t = ft+it*dt;
%           cost = cos(t);
%           sint = sin(t);
%           x1(c) = i1+a*cost*v1-b*sint*v2;
%           x2(c) = i2+a*cost*v2+b*sint*v1;
%         end
%         
%       end
% 
% 
% plot(x1,x2,'.y','MarkerSize',4.5)
% %         ylim([0 m2])
% axis ij
% axis equal

%
ev = interp.Sm.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));
    for i=1:mesh.num_param
        i1 = X(i);
        i2 = Y(i);
        if i2==Y(1)
    i2=Y(1+m1);
    end
        if X(i)==X(1)||X(i)==X(2)
        i1= X(3);
    else if X(i)==X(m1)||X(i)==X(m1-1)
            i1= X(m1-2);
        end
    end
        v = et.getEigenvectorV(i1,i2);
%         s = et.getEigenvalues(i1,i2);
s = ev(:,i2,i1);

        v1 = v(1); v2 = v(2);
        du = s(1); dv = s(2);
% if abs(s(1)-s(2))<1
% %     du=0.25; dv=1;
%     v1=1;v2=0;
% end
i1 = mod(i,m1);
if(i1==0) i1=m1-1;

else
    i1=i1-1; 
end
    i2 = ceil(i/m1)-1; 
    
    i1=mesh.param_x(i1+1);
          i2=YYY(i2+1);
        a = r*sqrt(dv); b = r*sqrt(du);
        for it=1:nt
          c=c+1;
          t = ft+it*dt;
          cost = cos(t);
          sint = sin(t);
          x1(c) = i1+a*cost*v1-b*sint*v2;
          x2(c) = i2+a*cost*v2+b*sint*v1;
        end
        
      end
 figure

plot(x1,x2,'.','MarkerSize',1)
        
axis ij
axis equal
xlim([-9 42])
ylim([0 75])      