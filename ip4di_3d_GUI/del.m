
tmp=log10(final.all_res_model(:,6));%remmber if is
tmp1=1000*(final.ip_model(:,6)-final.all_res_model(:,6))./(final.ip_model(:,6));


      min_x=min(final.param_x);
        max_x=max(final.param_x);

        min_y=min(final.param_y);
        max_y=max(final.param_y);

        min_z=min(final.param_z);
        max_z=max(final.param_z);
       xx=union(final.param_x, final.param_x); % take the unique x_line
        yy=union(final.param_y, final.param_y); %take the unique y_line
        zz=union(final.param_z, final.param_z); %take the unique y_line
        
        
        tmp_x=xx(2)-xx(1);% find the x step
        tmp_y=yy(2)-yy(1);% find the y_step
        tmp_z=zz(2)-zz(1);
        
        stepp=16;
xnew=[min_x:tmp_x/stepp:max_x];
ynew=[min_y:tmp_y/stepp:max_y];
znew=[0:tmp_z/stepp:max_z];



[xxx,yyy,zzz]=meshgrid(xnew,ynew,znew);
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


x_unique=unique(final.param_x);
y_unique=unique(final.param_y);
z_unique=unique(final.param_z);


%-----------------------Plots----------------------------------------------
    %ZI=griddata(param_x-neg_x,-param_y,tmp,XI,YI,'nearest');
figure    
subplot(2,1,1)
FF=TriScatteredInterp(final.param_x,final.param_y,final.param_z,tmp);
vi=FF(xxx,yyy,zzz);
slice(xxx,yyy,-zzz,vi,[],[],[-z_unique]);
shading flat
set(gca,'Box','On');
colorbar
colormap(map(1:17,:));

    title(['Amplitude. Iteration ',num2str(6),' RMS=> ',num2str(3.38),'%']);    
%     axis ('equal');
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    zlabel('Distance (m)');
%     xlim([min_x max_x]);
%     ylim([min_y max_y]);



  
subplot(2,1,2)
tmp2=0.7*(tmp1+60);
    FF=TriScatteredInterp(final.param_x,final.param_y,final.param_z,tmp2 );


vi=FF(xxx,yyy,zzz);
slice(xxx,yyy,-zzz,vi,[],[],[-z_unique]);
    colormap(map(1:17,:));
%     colormap(cmap);
 title(['Chargeability. Iteration ',num2str(6),' RMS=> ',num2str(3.34),'%']);    

shading ('flat');
set(gca,'Box','On');
 colorbar;

%     axis ('equal');
    xlabel('Distance (m)');
    ylabel('Distance (m)');
    zlabel('Distance (m)');


colormap(map(1:17,:))