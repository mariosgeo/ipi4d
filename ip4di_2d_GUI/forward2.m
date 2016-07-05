%% Elementary stuff
freq=50;
freq=2*freq;


xo=1;
zo=1;

xmax=50;
zmax=30;

dx=2;
dz=2;

[X Y] = meshgrid(xo:dx:xmax, zo:dz:zmax);


%%  Creadte SPEED model

xx=[xo:dx:xmax];
zz=[zo:dz:zmax];
xf=zeros(length(xx)*length(zz),1);
zf=zeros(length(xx)*length(zz),1);
SpeedImage=zeros(length(xx)*length(zz),1);
tmp1=1;
for i=1:length(xx)
    for j=1:length(zz)
        xf(tmp1)=xx(i);
        zf(tmp1)=zz(j);
        
        if zz(j)>=10 %&& j<=8 && i>=7 && i<=9
            SpeedImage(tmp1)=3;
        else
            SpeedImage(tmp1)=1;
        end
        
        
        tmp1=tmp1+1;
    end
end
SpeedImage=1000*SpeedImage; %m/s
VI=TriScatteredInterp(xf,zf,SpeedImage);


SpeedImage=VI(X,Y);

subplot(241); contourf(X,Y,SpeedImage); title('Real Model');colorbar;
% SpeedImage=1./SpeedImage;
%%
% Create input that has 4 columns, Sx, Sz, Rx, Rz
% Possible nodes that can act as Source-Receiver locations
% Create a p-shape source-receicer location

xgeo=unique(xf);
zgeo=unique(zf);
% left
tmp1=1;
for i=1:length(zgeo)
    g(tmp1,1)=min(xgeo);
    g(tmp1,2)=zgeo(i);
    tmp1=tmp1+1;
end
%surface
for i=1:length(xgeo);
    g(tmp1,1)=xgeo(i);
    g(tmp1,2)=min(zgeo);
    tmp1=tmp1+1;
end

% right
for i=1:length(zgeo)
    g(tmp1,1)=max(xgeo);
    g(tmp1,2)=zgeo(i);
    tmp1=tmp1+1;
end
g=union(g,g,'rows');


tmp11=1;
number_elec=length(g);
pos_elec=[1:1:number_elec];

for i=1:length(pos_elec)
    a=i;
    rest_electrodes=setdiff(pos_elec,a);
    s_size=size(rest_electrodes);
    for j=1:s_size(2)
        m=rest_electrodes(j);
        dd(tmp11,1)=a;
        dd(tmp11,2)=rest_electrodes(j);
        tmp11=tmp11+1;
    end
        
end


% assign coordinates to sources receivers
data=zeros   (      length(dd(:,1)),   4);
for i=1:length(dd(:,1))
    
        
    data(i,1)=g(dd(i,1),1); % Sx
    data(i,2)=g(dd(i,1),2); % Sz
    
    data(i,3)=g(dd(i,2),1); %Rx
    data(i,4)=g(dd(i,2),2); %Rz
end
    
    
num_mes=length(data(:,1));




%%
[AA,SLOW_in,XI,YI]=create_A_matrix(1./SpeedImage,X,Y,1,1)  ;   

[alltimes] = graphallshortestpaths(AA,'directed',false);    

B=ComputeBC(1,0); %petros function

new_x=unique(XI);
new_z=unique(YI);
x=reshape(XI,length(new_x)*length(new_z),1);
z=reshape(YI,length(new_x)*length(new_z),1);

J=zeros(num_mes,length(unique(X))*length(unique(Y)));
for i=1:num_mes
    SourcePoint=[data(i,1);data(i,2)]; % Actual Coordinates
    ReceiverPoint=[data(i,3);data(i,4)];


    % Find source ID
    tmpr1=find(x==SourcePoint(1,ite));
    tmpr2=find(z==SourcePoint(2,ite));
    source=intersect(tmpr1,tmpr2);
    % Find Recever ID
    tmpr1=find(x==ReceiverPoint(1,ite));
    tmpr2=find(z==ReceiverPoint(2,ite));
    receiver=intersect(tmpr1,tmpr2);

    % calculate path
    [dist,path,pred] = graphshortestpath(AA,source,receiver,'directed',false);
    path_x=zeros(length(path),1);
    path_z=zeros(length(path),1);
    fullpath=zeros(2*length(path),1);
    % new mesh grid
    tmp=1;
    for i=1:length(path)
        path_x(i)=x(path(i));
        fullpath(tmp)=path_x(i);
        tmp=tmp+1;
        path_z(i)=z(path(i));
        fullpath(tmp)=path_z(i);
        tmp=tmp+1;    
    end
    
    [S(i),T(i)]=pathlength(XI,YI,SLOW_in,fullpath,B);  % Petros function   
    temptimes=abs(alltimes(source,:)+alltimes(receiver,:)-alltimes(source,receiver));
    fresnel_times=zeros(1,length(new_x)*length(new_z));
    
    fresnel_time_ind=find(temptimes<=1/(freq)); % keep index of times we need
    fresnel_times(fresnel_time_ind)=1-freq*temptimes(fresnel_time_ind); %there are the weights
    
    W=TriScatteredInterp(x,z,fresnel_times');
    W=W(XI,YI);
    
    W_small=interp2(XI,YI,W,X,Y);
    J(i,:)=reshape(W_small,length(unique(X))*length(unique(Y)),1);
    J(i,:)=J(i,:).*(S(i)/sum(J(i,:))); 
    
    
end


final=[data T];

