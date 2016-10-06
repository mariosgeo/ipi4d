%
% Function AUTO_CONTOUR: plots the intermediate models reconstructed at each iteration
%
% Author: Marios Karaoulis, Colorado School of Mines

% Modified by Francois Lavoue', Colorado School of Mines, September 2015.
% Modifications include:
% 1) some explicit comments for more readability,
% 2) some change of notations. In particular:
%    - 'tmp'  is now 'amp' (amplitude)
%    - 'tmp1' is now 'phi' (phase)
% 3) simplification of the choice of the colormap
%    (now defined by user in 'plot_parameters.m' for
%     consistency with other plots throughout the code)

function y=auto_contour(log_flag,itr,ip_cnt,input,mesh,fem)

global big_part ip_part

figure
big_part=gca;
% cla (big_part);
figure
ip_part=gca;


if input.sip_flag==0 % Plots for DC and IP
    if input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==1) && input.sip_flag==0
        amp=log10(mesh.res_param1);%remmber if is
    elseif input.ip_flag==1 && ip_cnt==2 && input.sip_flag==0
        phi=1000*(mesh.res_param1-mesh.res_final_param)./(mesh.res_param1);   %FL: is this a phase?
        amp=log10(mesh.res_final_param);
    end

elseif input.sip_flag==1 % Plots for SIP
    amp=log10(abs(mesh.res_param1));
    phi=1000*atan2(imag(mesh.res_param1),real(mesh.res_param1));   
end


%-----------------------Plots----------------------------------------------

if input.plot_options.interp==0
% plot discretized model

   vx=unique(mesh.param_x); nx=length(vx);
   vy=unique(mesh.param_y); ny=length(vy);
   model=reshape(amp,[nx,nz]);
   model=transpose(model);
   imagesc(vx,vz,model);

elseif input.plot_options.interp==1
% plot interpolated model

   %ZI=griddata(param_x-neg_x,-param_y,amp,XI,YI,'nearest');
   ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y-mesh.param_ytopo,amp);
   ZII = ZI(mesh.xxx,mesh.yyy);
   contourf(big_part,mesh.xxx,mesh.yyy,ZII,17,'EdgeColor','none');
end


% tune figure
hleg=colormap(input.plot_options.cmap);
colorbar('peer',big_part);
title(big_part,['Amplitude. Iteration ',num2str(itr),' normalized RMS = ',num2str(fem.nrms)]);
axis (big_part,'equal');
xlabel(big_part,'Distance (m)');
ylabel(big_part,'Depth (m)');
xlim(big_part,[mesh.min_x mesh.max_x]);
ylim(big_part,[-mesh.max_y -mesh.min_y]);
  

% Plot the phase in case of IP or SIP data
if input.sip_flag==1 || input.ip_flag==1 && ip_cnt==2

    if input.plot_options.interp==0
    % plot discretized model

       vx=unique(mesh.param_x); nx=length(vx);
       vy=unique(mesh.param_y); ny=length(vy);
       model=reshape(phi,[nx,nz]);
       model=transpose(model);
       imagesc(vx,vz,model);

    elseif input.plot_options.interp==1
    % plot interpolated model

       %ZI=griddata(param_x-neg_x,-param_y,phi,XI,YI,'nearest');
       ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y-mesh.param_ytopo,phi);
       ZII = ZI(mesh.xxx,mesh.yyy);
       contourf(ip_part,mesh.xxx,mesh.yyy,ZII,17,'EdgeColor','none');
    end

    colormap(input.plot_options.cmap);

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
