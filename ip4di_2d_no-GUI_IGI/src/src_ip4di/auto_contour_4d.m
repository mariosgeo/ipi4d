function auto_contour_4d(itr,input,mesh,fem)

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

% loop over time-lapse datasets
for k=1:input.num_files
 figure()
 tmp2=tmp(:,k);    
 ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y,tmp2);
 ZII = ZI(mesh.xxx,mesh.yyy);
 contourf(mesh.xxx,mesh.yyy,ZII,'lineStyle','None');
 title(['Amplitude. Iteration ',num2str(itr),' Model ',num2str(k),' normalized RMS = ',num2str(fem.nrms)]);
 caxis ([min_c max_c]);
 colormap(cmap)
 xlabel('Distance (m)');
 ylabel('Depth (m)');
 xlim([mesh.min_x mesh.max_x]);
 ylim([-mesh.max_y -mesh.min_y]);
 
 if k==input.num_files
     colorbar('Location','EastOutside');
 end
 
end


if input.sip_flag==1 || input.ip_flag==1 && ip_cnt==2
   
   min_c=min(min(tmp1));
   max_c=max(max(tmp1));

   % loop over time-lapse datasets
   for k=1:input.num_files

      figure()
      tmp2=tmp1(:,k);    
      ZI=TriScatteredInterp(mesh.param_x,-mesh.param_y,tmp2);
      ZII = ZI(mesh.xxx,mesh.yyy);
      contourf(mesh.xxx,mesh.yyy,ZII,'lineStyle','None');
      title(['Phase. Iteration ',num2str(itr),' Model ',num2str(k)]);
      caxis ([min_c max_c]);
      colormap(cmap)
      xlabel('Distance (m)');
      ylabel('Depth (m)');
      xlim([mesh.min_x mesh.max_x]);
      ylim([-mesh.max_y -mesh.min_y]);
 
      if k==input.num_files
         colorbar('Location','EastOutside');
      end

   end

end

end   %end function auto_contour_4d
