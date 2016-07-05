function y=auto_contour(log_flag,itr,ip_cnt,input,mesh,fem)


global big_part ip_part
cla (big_part);




if input.sip_flag==0 % Plots for DC and IP
    if input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==1) &&input.sip_flag==0
        tmp=log10(mesh.res_param1);%remmber if is
    elseif input.ip_flag==1 && ip_cnt==2 &&input.sip_flag==0
        tmp1=1000*(mesh.res_param1-mesh.res_final_param)./(mesh.res_param1);
        tmp=log10(mesh.res_final_param);
    end
    
elseif input.sip_flag==1 % Plots for SIP
    tmp=log10(abs(mesh.res_param1));
    tmp1=1000*atan2(imag(mesh.res_param1),real(mesh.res_param1));   
end



%------------------Colormaps----------------------------------------------
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
    %ZI=griddata(param_x-neg_x,-param_y,tmp,XI,YI,'nearest');
    
    ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y-mesh.y_tmp,tmp);
    ZII = ZI(mesh.xxx,mesh.yyy);

if input.loke_colors==0    
    contourf(big_part,mesh.xxx,mesh.yyy,ZII,17,'EdgeColor','none');
    colormap(cmap)
else
    contourf(big_part,mesh.xxx,mesh.yyy,ZII,17,'EdgeColor','none');
    colormap(big_part,map(1:17,:));
end
    colorbar('peer',big_part);
    title(big_part,['Amplitude. Iteration ',num2str(itr),' RMS=> ',num2str(fem.rms_sum1)]);
    axis (big_part,'equal');
    xlabel(big_part,'Distance (m)');
    ylabel(big_part,'Depth (m)');
    xlim(big_part,[mesh.min_x mesh.max_x]);
    ylim(big_part,[-mesh.max_y -mesh.min_y]);
  

if input.sip_flag==1 || input.ip_flag==1 && ip_cnt==2
%     cla(ip_part);
    ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y-mesh.y_tmp,tmp1);
    ZII = ZI(mesh.xxx,mesh.yyy);

    contourf(ip_part,mesh.xxx,mesh.yyy,ZII,17,'EdgeColor','none');
%     colormap(map(1:17,:));
    colormap(cmap);
    if input.ip_flag==1
        title(ip_part,'Chargebility');
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
end