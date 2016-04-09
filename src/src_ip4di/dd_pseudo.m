function dd_pseudo(mesh,input)
global pp_plot

data=importdata(input.mes_in);
input.ax=data(:,1);
input.bx=data(:,3);
input.mx=data(:,5);
input.nx=data(:,7);



rr=data(:,9);
mag=rr;
if input.sip_flag==1
    ii=data(:,10);
    array_model_data=complex(rr,ii);
    
    mag=abs(array_model_data);
    phi=1000*atan(imag(array_model_data)./real(array_model_data));
end



for i=1:input.num_mes
    
    % find A-B
   x_1=abs(input.ax(i)+input.bx(i))/2;
   x_2=abs(input.mx(i)+input.nx(i))/2;
  
   x_center(i)=abs(x_1+x_2)/2;
   y_center(i)=abs(x_1-x_2);
   
    
end

min_x=min( x_center);
max_x=max( x_center);

min_y=min(y_center);
max_y=max(y_center);
xx=union( x_center,  x_center); % take the unique x_line
yy=union(y_center,y_center); %take the unique y_line
tmp_x=xx(2)-xx(1);% find the x step
tmp_y=yy(2)-yy(1);% find the y_step

xnew=[min_x:tmp_x/16:max_x];
ynew=[0:tmp_y/16:max_y];

ynew=ynew; %Anapoda giati einai arnitika
[xxx,yyy]=meshgrid(xnew,ynew);

[XI YI]=meshgrid(xx,yy);

ZI=TriScatteredInterp(x_center',y_center',(mag));
ZII = ZI(xxx,yyy);



contourf(pp_plot,xxx,-yyy,ZII);
shading interp
title(pp_plot,'Amplitude Ohm.m');
colorbar('peer',pp_plot);

if input.sip_flag==1
figure

ZI=TriScatteredInterp(x_center',y_center',(phi));
ZII = ZI(xxx,yyy);



contourf(pp_plot,xxx,-yyy,ZII);
shading interp
title('Phase mrad');
shading interp
colorbar
end
end