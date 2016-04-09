function preprocessing

%global mes_in
global act

input.num_files

cmap=cptcmap('GMT_seis');
cmap=cmap(end:-1:1,:); %MARIOS ADD IN


tmp=importdata(mes_in);
num_files=length(tmp);

% just useles...
%[misc,tin,misc,misc,misc,tmp15,misc,nitr,misc,tmp16_lag,misc,tpar,misc,tmes]=...
                                    textread(tmp{1},'%s %d %s %s %s %f %s %d %s %f %s %d %s %d',1);

  num_param=tpar;                              
                                
try
    for i=1:num_files
        data_i(:,:,i)=importdata(tmp{i},' ',1);
    end
catch
    disp('this version does not support files with different number of measurments on each file. Future wil...');
    return
end
    tmp=data_i(:,1,1);
    x=tmp.data(1:num_param,1,1);
    y=tmp.data(1:num_param,2,1);
    
    num_param=length(x);
    
    res_model(:,1)=tmp.data(1:num_param,3,1);
    
    sip_flag=0;
    try %if in first file, there is 10th column, then check for IP or SIP
    imag_model(:,1)=tmp.data(1:num_param,4,1);
    sip_flag=1;
    end
    
    
    for i=1:num_files
        tmp=data_i(:,1,i);
        res_model(:,i)=tmp.data(1:num_param,3);
       if sip_flag==1; imag_model(:,i)=tmp.data(1:num_param,4);end
    end
        
    

[xi,yi]=meshgrid(unique(x),unique(y));
l_final=[];
% compeare two files per time
% figure
for i=1:num_files-1

    % First search for resistivity (or amplitude)
    tact(:,1)=abs(log10(res_model(:,i+1))-log10(res_model(:,i)))+eps;
%     tact(:,1)=res_model(:,i+1)./res_model(:,i);
    SP_max=max(tact(:,1));
    SP_min=min(tact(:,1));
    L1=zeros(num_param,1);

    for j=1:num_param
        L1(j)=log10(lagrn_min) +  ( (log10(lagrn_max) - log10(lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(tact(j,1)) - log10(SP_min));    
        L1(j)=10^(L1(j));
    end

    L1=lagrn_max-L1+lagrn_min;

    subplot(num_files,sip_flag+1,i);
 tact(:,1)=res_model(:,i+1)-res_model(:,i);
%     FF=TriScatteredInterp(x,y,L1);
    FF=TriScatteredInterp(x,y,tact(:,1));
    vi=FF(xi,yi);

    contourf(xi,yi,vi,17,'EdgeColor','None');
    colormap(cmap);
    colorbar
    axis equal
    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);
    caxis([-100 100]);
    % If we have phase
    if sip_flag==1
       tact(:,2)=abs(imag_model(:,i+1)-imag_model(:,i));
       SP_max=max(tact(:,2));
       SP_min=min(tact(:,2));
       L2=zeros(num_param,1);
        for j=1:num_param
            L2(j)=log10(lagrn_min) +  ( (log10(lagrn_max) - log10(lagrn_min) ) / (log10(SP_max) - log10(SP_min)) )*(log10(tact(j,2)) - log10(SP_min));    
            L2(j)=10^(L2(j));
        end       
       L2=lagrn_max-L2+lagrn_min;
       subplot(sip_flag+1,num_files,i+num_files);
       FF=TriScatteredInterp(x,y,L2);
       vi=FF(xi,yi);

       contourf(xi,yi,vi,17,'EdgeColor','None');
       colormap(cmap);
       colorbar 
%        axis equal
       xlim([min(x) max(x)]);
       ylim([min(y) max(y)]);
    end
    
    % Keep final 
    LL=zeros(num_param,1);
    if sip_flag==1
        for j=1:length(L1);
            LL(j)=min(L1(j),L2(j));
        end
    end
    
    
    if sip_flag==1
        l_final=[l_final;LL];
    else
        l_final=[l_final;L1];
    end
    
end
    
 act=zeros(num_files*num_param,   num_files*num_param);
 for i=1:length(l_final)
    act(i,i)=l_final(i);
 end


