function auto_contour_4d(itr,input,mesh,fem)

global tim
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



if input.sip_flag==0 % Plots for DC and IP
    if input.ip_flag==0 || (input.ip_flag==1 && ip_cnt==1) &&input.sip_flag==0
        tmp=log10(mesh.d4_res_param1);%remmber if is
    elseif input.ip_flag==1 && ip_cnt==2 &&input.sip_flag==0
% %         tmp=1000*(mesh.d4.res_param2-mesh.res_final_param)./(mesh.res_param2);
% %         tmp1=log10(mesh.res_final_param);
    end
    
elseif input.sip_flag==1 % Plots for SIP
    tmp=log10(abs(mesh.d4_res_param1));
    tmp1=1000*atan2(imag(mesh.d4_res_param1),real(mesh.d4_res_param1));   
end





min_c=min(min(tmp));
max_c=max(max(tmp));
for k=1:input.num_files



 tmp2=tmp(:,k);    
 ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y,tmp2);
 ZII = ZI(mesh.xxx,mesh.yyy);
 contourf(tim(k),mesh.xxx,mesh.yyy,ZII);
 title(tim(k),['Amplitude. Iteration ',num2str(itr),' Model ',num2str(k),' RMS',num2str(fem.rms_sum1)]);
 caxis (tim(k),[min_c max_c]);
 colormap(cmap)
 xlabel(tim(k),'Distance (m)');
 ylabel(tim(k),'Depth (m)');
 xlim(tim(k),[mesh.min_x mesh.max_x]);
 ylim(tim(k),[-mesh.max_y -mesh.min_y]);
 
 if k==input.num_files
     colorbar('Location','EastOutside');
 end
 
end




if input.sip_flag==1 || input.ip_flag==1 && ip_cnt==2
   
    min_c=min(min(tmp1));
max_c=max(max(tmp1));
    for k=1:input.num_files



 tmp2=tmp1(:,k);    
 ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y,tmp2);
 ZII = ZI(mesh.xxx,mesh.yyy);
 contourf(tim(k+input.num_files),mesh.xxx,mesh.yyy,ZII);
 title(tim(k+input.num_files),['Phase. Iteration ',num2str(itr),' Model ',num2str(k)]);
 caxis (tim(k+input.num_files),[min_c max_c]);
 colormap(cmap)
 xlabel(tim(k+input.num_files),'Distance (m)');
 ylabel(tim(k+input.num_files),'Depth (m)');
 xlim(tim(k+input.num_files),[mesh.min_x mesh.max_x]);
 ylim(tim(k+input.num_files),[-mesh.max_y -mesh.min_y]);
 
 if k==input.num_files
     colorbar('Location','EastOutside');
 end
 
    end

end

end