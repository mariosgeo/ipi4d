%
% SUBROUTINE TO PLOT MODEL ON INPUT GRID
%

function plot_model_on_input_grid(input,mesh,model)

 figure
 hold on

 % reshape as rectangular matrix
 model=reshape(model,[input.n1,input.n2]);

 % input vectors
 vx1=[0:input.n2-1]*input.hstep;    %row
 vz1=[0:input.n1-1]'*input.hstep;   %column

 if input.plot_options.cmplx_flag==1
    imagesc(vx1,vz1,real(model))
    hleg=colorbar;
    ylabel(hleg,'Resistivity (\Omega.m), real part')

 elseif input.plot_options.cmplx_flag==2
    imagesc(vx1,vz1,imag(model))
    hleg=colorbar;
    ylabel(hleg,'Resistivity (\Omega.m), imag. part')

 elseif input.plot_options.cmplx_flag==3
    imagesc(vx1,vz1,abs(model));
    hleg=colorbar;
    ylabel(hleg,'Resistivity (\Omega.m), amplitude')
    input.plot_options.caxis=input.plot_options.caxis_amp;
    input.plot_options.axis_tics=input.plot_options.axis_tics_amp;

 elseif input.plot_options.cmplx_flag==4
    imagesc(vx1,vz1,atan2(imag(model),real(model))*1000);
    hleg=colorbar;
    ylabel(hleg,'Resistivity, phase (mrad)')
    input.plot_options.caxis=input.plot_options.caxis_phi;
    input.plot_options.axis_tics=input.plot_options.axis_tics_phi;
 end

 plot(mesh.orig_probe_x,-mesh.orig_probe_z,'o');
 set(gca,'YDir','reverse');
 axis equal;

 title(input.plot_options.label_title)
 xlabel('x (m)')
 ylabel('z (m)')

 xlim([0,(input.n2-1)*input.hstep])
 ylim([0,(input.n1-1)*input.hstep])

if(numel(input.plot_options.caxis)>0); caxis(input.plot_options.caxis); end
if(numel(input.plot_options.axis_tics)>0); set(hleg,'YTick',input.plot_options.axis_tics); end

 % reshape as single vector
 model=model(:);

 % clean
 clear vx1 vz1
