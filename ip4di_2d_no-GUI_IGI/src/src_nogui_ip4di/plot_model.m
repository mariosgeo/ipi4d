%
% Function PLOT_MODEL: simply plot model when mesh is not known
%
% Derived from function auto_contour.m by M. Karaoulis.
%
% Author: Francois Lavoue', Colorado School of Mines
% Version: October 15, 2015.

function plot_model(vx,vz,model,input)

figure

if input.ip_flag==0 && input.sip_flag==0 % Plots for DC
   amp=log10(model);

elseif input.ip_flag==1
   disp('FUNCTION plot_model.m IS NOT VALID FOR IP DATA.')
   return;   %get out of function

elseif input.sip_flag==1 % Plots for SIP

    amp=log10(abs(model));
    phi=1000*atan2(imag(model),real(model));   
end


%------------------ Res2Dinv colormap ----------------------------------------------
map(1,:)=[0 0 128/255];
map(2,:)=[0 0 170/255];
map(3,:)=[0 0 211/255];
map(4,:)=[0 0 255/255];
map(5,:)=[0 128/255 255/255];
map(6,:)=[0 255/255 255/255];
map(7,:)=[0 192/255 128/255];
map(8,:)=[0 255/255 0];
map(9,:)=[0 128/255 0];
map(10,:)=[128/255 192/255 0];
map(11,:)=[255/255 255/255 0];
map(12,:)=[191/255 128/255 0];
map(13,:)=[255/255 128/255 0];
map(14,:)=[255/255 0 0];
map(15,:)=[211/255 0 0];
map(16,:)=[132/255 0 64/255];
map(17,:)=[96/255 0/255 96/255];
map(18,:)=[255/255 255/255 255/255];


cmap=cptcmap('GMT_seis');
cmap=cmap(end:-1:1,:); %MARIOS ADD IN


%-----------------------Plots----------------------------------------------

% mesh
[MX,MZ]=meshgrid(vx,vz);
param_x=MX(:);
param_z=MZ(:);

% INTERPOLATE MODEL
%ZI=griddata(param_x-neg_x,-param_y,amp,XI,YI,'nearest');
ZI=TriScatteredInterp(param_x,param_z,amp);
ZII = ZI(MX,MZ);

if input.loke_colors==0    
%plot with GMT_seis colormap
    contourf(big_part,MX,MZ,ZII,17,'EdgeColor','none');
    colormap(cmap)
else
%plot with Res2Dinv colormap
    contourf(big_part,MX,MZ,ZII,17,'EdgeColor','none');
    colormap(big_part,map(1:17,:));
end


% Plot the amplitude
colorbar('peer',big_part);
title(big_part,['Amplitude. Iteration ',num2str(itr),' normalized RMS = ',num2str(fem.nrms)]);
axis (big_part,'equal');
xlabel(big_part,'Distance (m)');
ylabel(big_part,'Depth (m)');
xlim(big_part,[x_min x_max]);
ylim(big_part,[z_min z_max]);
  

% Plot the phase in case of IP or SIP data
if input.sip_flag==1 || input.ip_flag==1 && ip_cnt==2
%     cla(ip_part);
    ZI=TriScatteredInterp(param_x,param_z,phi);
    ZII = ZI(MX,MZ);

    contourf(ip_part,MX,MZ,ZII,17,'EdgeColor','none');
%     colormap(map(1:17,:));
    colormap(cmap);
    if input.ip_flag==1
       title(ip_part,'Chargeability');
    elseif input.sip_flag==1
       title(ip_part,'Phase');
    end
    axis (ip_part,'equal');
    xlabel(ip_part,'Distance (m)');
    ylabel(ip_part,'Depth (m)');
    xlim(ip_part,[mesh.min_x mesh.max_x]);
    ylim(ip_part,[-mesh.max_y -mesh.min_y]);
    colorbar('peer',ip_part);
   % drawnow;
end

end   %end function
