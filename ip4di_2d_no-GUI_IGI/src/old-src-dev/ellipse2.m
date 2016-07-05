%
% PLOT ELLIPSES REPRESENTING TENSORS
%

%% TUNING PARAMETERS
r = 0.01;   %parameter to scale the size of ellipse                         
nt = 101;   %discretize theta for plotting
usamp=1;    %under-sample grid (don't plot ellipses at all grid points
            %  to improve computation time and readability)

marker_size=10
%%END TUNING PARAMETERS


dt = 2.0*pi/(nt-1);
ft = 0;

x1=zeros(nt,1);
x2=zeros(nt,1);

c=0;
XXX=unique(mesh.param_x);
YYY=unique(mesh.param_y);

%% RAPPEL CONVENTIONS:
% X = vector of size num_param with all INDICES of grid cells in the x-direction
% Y =   "    "   "      "       "    "     "    "   "     "   "   "  y-direction
% mesh.param_x = vector of size num_param with all LOCATIONS of grid cells in the x-direction
% mesh.param_y =   "    "   "      "       "    "      "     "   "     "   "   "  y-direction
% XXX = vector of size M1 with LOCATIONS of cells in the x-direction
% YYY =   "    "   "   M2  "       "     "    "   "   "  y-direction


%FL DEBUG
val_n1=n1
val_n2=n2
val_m1=m1
val_m2=m2

%
ev = interp.Sm.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));
%ev = interp.Sm.getEigenvalues(et,zeros(m2,m1),zeros(m2,m1));
%FL: ev is over-sized if under-sampling is used...

%for i=1:mesh.num_param
%    i1 = X(i);
%    i2 = Y(i);

for i2=1:usamp:m2

  for i1=1:usamp:m1

    if i2==Y(1)
       i2=Y(1+m1);
    end

    %if X(i)==X(1)||X(i)==X(2)
    if i1==X(1) || i1==X(2)
       i1= X(3);
    else 
       %if X(i)==X(m1)||X(i)==X(m1-1)
       if i1==X(m1) || i1==X(m1-1)
          i1= X(m1-2);
       end
    end

    if i1<m1 && i2<m2
       %s = et.getEigenvalues(i1,i2);
       v = et.getEigenvectorV(i1,i2);
       s = ev(:,i2,i1);

       v1 = v(1);
       v2 = v(2);
       du = s(1);
       dv = s(2);

% if abs(s(1)-s(2))<1
% %     du=0.25; dv=1;
%     v1=1;v2=0;
% end

       %i1 = mod(i,m1);
       %i2 = ceil(i/m1)-1;
       if (i1==0)
          i1=m1-1;
       else
          i1=i1-1; 
       end
    
       vi1=XXX(i1+1);
       vi2=YYY(i2+1);
       a = r*sqrt(dv);
       b = r*sqrt(du);

       for it=1:nt
           c=c+1;
           t = ft+it*dt;
           cost = cos(t);
           sint = sin(t);
           x1(c) = vi1+a*cost*v1-b*sint*v2;
           x2(c) = vi2+a*cost*v2+b*sint*v1;
       end

    end   %end if i1<m1 && i2<m2
  end   %end i1

end   %end i2
%end   %for i=1:mesh.num_param


% plot
figure
plot(x1,x2,'.','MarkerSize',marker_size)
        
axis ij
axis equal
xlim([0 max(mesh.param_x)])
ylim([0 max(mesh.param_y)])

% FL DEBUG
disp('Nb of ellipses plotted = ')
disp(c/nt)
