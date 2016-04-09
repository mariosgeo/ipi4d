%
% PLOT ELLIPSES REPRESENTING TENSORS
%

disp(' ')
disp('================================================')
disp('=              ENTER ELLIPSE.m                 =')
disp('= Plot ellipses representing structure tensors =')
disp('================================================')
disp(' ')

%% TUNING PARAMETERS

%parameter to scale the size of ellipse
%r = 0.0075;   %if using the true sandbox structure
%r = 0.05;     %if using the gradient image
r = 0.01;     %if using the gradient image and p0=0

nt = 101;   %discretize theta for plotting
usamp=1;    %under-sample grid (don't plot ellipses at all grid points
            %  to improve computation time and readability)

% plot options
marker_size=10;

%%END TUNING PARAMETERS


% parameter-variables for the ellipses
dt = 2.0*pi/(nt-1);
ft = 0;

% init. ellipse points locations
x1=zeros(nt,1);
x2=zeros(nt,1);

% init. count of points
c=0;



% LOOP OVER ROWS IN MODEL (y-indices)
for i2=1:usamp:m2

   % LOOP OVER COLUMNS IN MODEL (x-indices)
   for i1=1:usamp:m1

       % compute indices (k1,k2) in TI that corresponds to inversion cell (i1,i2)
       k1=ceil(XXX(i1)/input.hstep)+1;
       k2=ceil(YYY(i2)/input.hstep)+1;

       % Get eigenvalues of structure tensors
       %%(from ev computed in initcm2.m with JZ's interp.Sm.getEigenvalues)
       s = ev(:,k2,k1);

       % Get eigenvalues of structure tensors
       % (using Dave Hale's routine, gives the same results)
       %s = et.getEigenvalues(k1-1,k2-1);

       %FL: do NOT set the 2nd eigenvalue to 1, as done befoe in initcm2,
       %    because it leads to anisotropic smoothing even in isotropic regions.
       %s(2)=1;

       % Get the eigenvector V (parallel to structures) at point (k1,k2) in TI
       % (be careful to Java indices, starting at 0)
       v = et.getEigenvectorV(k1-1,k2-1);

       v1 = v(1);
       v2 = v(2);
       du = s(1);
       dv = s(2);

       % location of ellipse center    
       vi1=XXX(i1);
       vi2=YYY(i2);

       % sizes of ellipse axis
       a = r*sqrt(dv);   %FL2JZ: why sqrt?
       b = r*sqrt(du);

       for it=1:nt
           c=c+1;
           t = ft+it*dt;
           cost = cos(t);
           sint = sin(t);

           % parametric equations of the ellipse,
           % with large axis a = sqrt(av)*V parallel to structures
           %  and small axis b = sqrt(au)*U perpendicular to structures
           x1(c) = vi1 + a*cost*v1 - b*sint*v2;
           x2(c) = vi2 + a*cost*v2 + b*sint*v1;
       end

   end   %end i1

end   %end i2


% plot ellipses
if input.plot_ellipses==1
   figure
   plot(x1,x2,'.','MarkerSize',marker_size)
        
   axis ij
   axis equal
   xlim([0 max(mesh.param_x)])
   ylim([0 max(mesh.param_y)])
   hold on

else
% save ellipse coordinates for later plot

   Mts=[x1 x2];
   save('ellipses.dat','Mts','-ascii')

end

