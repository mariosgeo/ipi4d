 function [ctc,L1]=acb3d(input,mesh,fem)

global big_part ip_part
L1=zeros(mesh.num_param,1);
SP=zeros(mesh.num_param,1);

% tmp_lag=lagrn_min*mesh.jacscale;
   JTJ=fem.array_jacobian.'*fem.array_jacobian;

   %ctc=c'*c;
   resolution_matrix=(JTJ + 0.005*mesh.ctc)\(fem.array_jacobian.'*fem.array_jacobian);



   for i=1:mesh.num_param
       cur1=i;
       tmp999=0.0;
       for j=1:mesh.num_param
           cur2=j;
           if (cur1>mesh.num_param) cur1=i-mesh.num_param; end
           if (cur2>mesh.num_param) cur2=j-mesh.num_param; end
           
           	x1_dis=mesh.param_i(cur1,1);
    		x2_dis=mesh.param_i(cur2,1);
		
        	y1_dis=mesh.param_i(cur1,2);
            y2_dis=mesh.param_i(cur2,2);
            
            z1_dis=mesh.param_i(cur1,3);
            z2_dis=mesh.param_i(cur2,3);
            
                        
	    	weight=( (x1_dis - x2_dis)*(x1_dis - x2_dis) + (y1_dis - y2_dis)*(y1_dis - y2_dis) +(z1_dis - z2_dis)*(z1_dis - z2_dis)); %/*Not sqrt because in SP i need weiht*weight*/
           	
            tmp999=tmp999+ weight*( (1-mesh.S(i,j))*resolution_matrix(i,j) ) * ( (1-mesh.S(i,j))*resolution_matrix(i,j) ); 
       end
       
       SP(i)=tmp999;
   end

SP_max=max(SP);
SP_min=min(SP);


for i=1:mesh.num_param
	L1(i)=log10(input.lagrn_min) +  ( (log10(input.lagrn_max) - log10(input.lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(SP(i)) - log10(SP_min));    
    L1(i)=10^(L1(i));
end

minimum_security=0.5;
SP_security=minimum_security*( tmp999/mesh.num_param);


for i=1:mesh.num_param
    tmp999=0;
    for j=1:mesh.num_param
        tmp999=tmp999+resolution_matrix(i,j);
    end
    if (tmp999 <minimum_security)
        SP(i)=SP_max;
    end
end
for i=1:mesh.num_param
    if SP(i)>1
        SP(i)=1;
    end
end


for i=1:mesh.num_param
	L1(i)=log10(input.lagrn_min) +  ( (log10(input.lagrn_max) - log10(input.lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(SP(i)) - log10(SP_min));    
    L1(i)=10^(L1(i));
end  



LL=zeros(mesh.num_param,mesh.num_param);
for i=1:mesh.num_param
    LL(i,i)=L1(i);
end
% mesh.ctc=mesh.c'*LL*mesh.c;
% ctc=mesh.cx'*LL*mesh.cx + mesh.cy'*LL*mesh.cy + mesh.cz'*LL*mesh.cz;


ctc=zeros(mesh.num_param,mesh.num_param);
%update now ctc for export
for i=1:mesh.num_param
    for j=1:mesh.num_param
        ctc(i,j)=mesh.ctc(i,j)*real(L1(i)); %real in cases of complex numbers
    end
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


x_unique=unique(mesh.param_x);
y_unique=unique(mesh.param_y);
z_unique=unique(mesh.param_z);

min_x=min(x_unique);
max_x=max(x_unique);

min_y=min(y_unique);
max_y=max(y_unique);

min_z=min(z_unique);
max_z=max(z_unique);

[xi yi zi]=meshgrid(x_unique,y_unique,z_unique);
abs_up=max(log10(abs(mesh.res_param1)));
abs_down=min(log10(abs(mesh.res_param1)));
%-----------------------Plots----------------------------------------------
    %ZI=griddata(param_x-neg_x,-param_y,tmp,XI,YI,'nearest');
    
FF=TriScatteredInterp(mesh.param_x,mesh.param_y,mesh.param_z,(real(L1)));
vi=FF(xi,yi,zi);
slice(ip_part,xi,yi,-zi,vi,[],[],[-z_unique]);


 end