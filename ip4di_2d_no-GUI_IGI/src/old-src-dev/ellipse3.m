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
r = 0.008;   %parameter to scale the size of ellipse                         
nt = 101;   %discretize theta for plotting
usamp=1;    %under-sample grid (don't plot ellipses at all grid points
            %  to improve computation time and readability)

% plot
marker_size=10;
%%END TUNING PARAMETERS


dt = 2.0*pi/(nt-1);
ft = 0;

x1=zeros(nt,1);
x2=zeros(nt,1);

c=0;

% already defined in initcm2.m...
XXX=unique(mesh.param_x);
YYY=unique(mesh.param_y);

%% RAPPEL CONVENTIONS:
% X = vector of size num_param with all INDICES of grid cells in the x-direction
% Y =   "    "   "      "       "    "     "    "   "     "   "   "  y-direction
% mesh.param_x = vector of size num_param with all LOCATIONS of grid cells in the x-direction
% mesh.param_y =   "    "   "      "       "    "      "     "   "     "   "   "  y-direction
% XXX = vector of size M1 with LOCATIONS of cells in the x-direction
% YYY =   "    "   "   M2  "       "     "    "   "   "  y-direction


%
%ev = interp.Sm.getEigenvalues(et,zeros(m2,m1),zeros(m2,m1));
ev = interp.Sm.getEigenvalues(et,zeros(n2,n1),zeros(n2,n1));


%for i=1:mesh.num_param
%    i1 = X(i);
%    i2 = Y(i);


% LOOP OVER ROWS
for i2=1:usamp:m2

  % skip first row
  if i2==Y(1)
     i2=Y(1+m1);
  end

%FL DEBUG
val_i2=i2

  for i1=1:usamp:m1

%FL DEBUG
val_i1=i1

    %if X(i)==X(1)||X(i)==X(2)
    if i1==X(1) || i1==X(2)
    % skip first two columns
       i1= X(3);
    else 
    % skip last two columns
       %if X(i)==X(m1)||X(i)==X(m1-1)
       if i1==X(m1) || i1==X(m1-1)
          i1= X(m1-2);
       end
    end

%FL DEBUG
val_i1=i1

    %ensure we don't plot any ellipse at the edges
    if i1<m1 && i2<m2

%FL DEBUG
val_i2=i2
ii2=ceil(YYY(i2)/input.hstep)+1
val_YYY_i2=YYY(i2)
val_YY_ii2=(ii2-1)*input.hstep
disp(' ')

val_i1=i1
ii1=ceil(XXX(i1)/input.hstep)+1
val_XXX_i1=YYY(i1)
val_XX_ii1=(ii1-1)*input.hstep

if val_YYY_i2==val_YY_ii2 && val_XXX_i1==val_XX_ii1
   disp('------------')
else
   disp('ERROR')
%   break
end

       %s = et.getEigenvalues(i1,i2);
       v = et.getEigenvectorV(ii1,ii2);
       s = ev(:,ii2,ii1);

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


% plot ellipses
if input.plot_ellipses==1
   figure
   plot(x1,x2,'.','MarkerSize',marker_size)
        
   axis ij
   axis equal
   xlim([0 max(mesh.param_x)])
   ylim([0 max(mesh.param_y)])
   hold on

   % FL DEBUG
   disp('Nb of ellipses plotted = ')
   disp(c/nt)

else
% save ellipse coordinates for later plot

   Mts=[x1 x2];
   save('ellipses.dat','Mts','-ascii')

end

