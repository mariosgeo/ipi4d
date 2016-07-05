%
% Function PLOT_JACOBIAN: plot Jacobian matrix
%                         as a movie wrt to electrode combinations.
%
% Derived from MK's function pushbutton3_Callback
% in IP4DI 'forward_modelling.m'.
%
% Author: Francois Lavoue, Colorado School of Mines
% Version: October 15, 2015

function plot_Jacobian(input,mesh,fem)

figure

for l=1:input.num_mes
   cla
   %set(handles.edit1,'String',num2str(l));

   F = TriScatteredInterp(mesh.param_x,-mesh.param_y,fem.array_jacobian(l,:)');
   qz=F(mesh.xxx,mesh.yyy);
   %contourf(xxx,yyy,qz);
   contourf(mesh.xxx,mesh.yyy,qz);
   %shading interp   % commented FL

   colorbar
   axis equal;% off; 

   hold on
   plot(input.ax(l),-input.az(l),'ko','LineWidth',4);

   text(input.ax(l),-input.az(l)+mesh.probe_spacing,'A','FontSize',22);

   plot(input.bx(l),-input.bz(l),'ko','LineWidth',4);
   text(input.bx(l),-input.bz(l)+mesh.probe_spacing,'B','FontSize',22);

   plot(input.mx(l),-input.mz(l),'ko','LineWidth',4);
   text(input.mx(l),-input.mz(l)+mesh.probe_spacing,'M','FontSize',22);

   plot(input.nx(l),-input.nz(l),'ko','LineWidth',4);
   text(input.nx(l),-input.nz(l)+mesh.probe_spacing,'N','FontSize',22);

   title('JACOBIAN MATRIX')
   pause(0.01)
end
