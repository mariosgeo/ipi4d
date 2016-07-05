%
% PLOT_RESOLUTION: compute and plot resolution matrix (with data contribution only)
%
% Derived from function Resolution_button_Callback
% --- Executes on button press in Resolution_button.
% in Marios Karaoulis' forward_modelling.m.
% F. Lavoue', September 5, 2015

function plot_resolution(input,mesh,fem)

% Dear user since I have already calculated the jacobian (sensitivity)
% it's now easy for you to calculate the resolution matrix. Since this
% program is meant just for forward modelling, I didn't bother calculating
% a decent smoothness matrix. I just use DAMPED LEAST SQUARES. You can
% modify the equation below, so you can insert your smoothness matrix.
% In case someone needs a smoothness matrix, please look at the accompaning
% m files below and uncomment below lines, and and change the
% eye(num_param) with ctc (assuming c is the name of the smoothness).

% global num_of_bor
%if num_of_bor==1 || num_of_bor==2
%    smooth_mtx_horizontal;
%    
%elseif num_of_bor==0
%    smooth_mtx_surface;   
%end

% R = (J'J+eps*I)^-1 * J'J
resolution_matrix=(fem.array_jacobian.'*fem.array_jacobian + 0.005*eye(mesh.num_param))\(fem.array_jacobian.'*fem.array_jacobian);

% I chose to plot on the diagonal of the matrix
for i=1:mesh.num_param
%    resolution_matrix2(i,1)=resolution_matrix(i,i);
    resolution_matrix2(i,1)=log10(resolution_matrix(i,i));
end


% And now plot it...
figure
F = TriScatteredInterp(mesh.param_x,mesh.param_y,resolution_matrix2);
qz=F(mesh.xxx,-mesh.yyy);
contourf(mesh.xxx,-mesh.yyy,qz);
%shading interp

title('RESOLUTION MATRIX')

axis equal;
set(gca,'YDir','reverse');

xlim([0,(input.n2-1)*input.hstep])
ylim([0,(input.n1-1)*input.hstep])

hleg=colorbar;
%set(hleg,'YTick',[0,0.5,0.75,1])
set(hleg,'YTick',[-3,-2,-1,0])
set(hleg,'YTickLabel',[0.001,0.01,0.1,1])

