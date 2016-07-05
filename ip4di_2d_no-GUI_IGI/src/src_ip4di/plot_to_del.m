cmap=cptcmap('GMT_seis');
cmap=cmap(end:-1:1,:); %MARIOS ADD IN
x=final.param_x;
y=final.param_y;

num_param=length(x);

[xi,yi]=meshgrid(unique(x),unique(y));
res_model=final.d4_res_param1(:,:,5);
for i=1:2
    
    % First search for resistivity (or amplitude)
%     tact(:,1)=abs(log10(res_model(:,i+1))-log10(res_model(:,i)));
    tact(:,1)=res_model(:,i+1)./res_model(:,i);
    
    subplot(1,2,i);
       FF=TriScatteredInterp(x,y,tact(:,1));
    vi=FF(xi,yi);

    contourf(xi,-yi,vi,17,'EdgeColor','None');
    colormap(cmap);
    colorbar
%     axis equal
    xlim([min(x) max(x)]);
    ylim([min(-y) max(-y)]);
    xlabel('Distance (m)');
    ylabel('Depth (m)');
    title('T2/T1');
    caxis([0.85 2.4]);
    
end
    