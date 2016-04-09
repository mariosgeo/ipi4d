model2=inline('p(4)-0.5*p(1).*(1-sinh(p(2).*log(x*p(3)))./(cosh(p(2).*log(x*p(3)))+sin(pi/2*(1-p(2)))))','p','x');

for i=1:M
    Y=imag(final_model(i,:));
    Y2=real(final_model(i,:));
    [beta2(i,:),r2,J2] =nlinfit2(w',Y2',model2,[max(Y2),real(beta(i,:))]);   
end

%% petrophysical parameters for the 4 blocks of the synthetic model
sigma_inf=[ 17.7200  116.1100   69.6700   26.0500]*1e-4;
sigma_zero=[17.5100   95.1100   57.0700   25.8400]*1e-4;
Mn=[0.2100   21.0000   12.6000    0.2100]*1e-4;
tau0=[3 0.2 0.02 25];
c=[0.5 0.5 0.25 0.45];
trueF=[6 9 15 4];

%% simulated frequencies
f=[0.0010    0.0032    0.0100    0.0316    0.1000    0.3162    1.0000   3.1623   10.0000   31.6228  100.0000];
w=2*pi*f;

%% apply petrophysical relationships to define sigma_r and sigma_i for each frequency
for i =1:length(w)
    sigma_r(:,i)=sigma_inf-0.5*Mn.*(1-sinh(c.*log(w(i)*tau0))./(cosh(c.*log(w(i)*tau0))+sin(pi/2*(1-c))));
    sigma_i(:,i)=-0.5*Mn.*cos(pi/2*(1-c))./(cosh(c.*log(w(i)*tau0))+sin(pi/2*(1-c)));
    amplitude(:,i)=1./((sigma_r(:,i).^2+sigma_i(:,i).^2).^0.5);
    phase(:,i)=atan(-sigma_i(:,i)./sigma_r(:,i))*1000;   % in mrad
%   phase2(:,i)=atan2(-sigma_i(:,i),sigma_r(:,i))*1000;
end
%%

%% build synthetic model for each frequency
for i=1:length(f)
    model4_chsip(find(model4(:,3)==1000),3)=amplitude(4,i);
    model4_chsip(find(model4(:,3)==1000),4)=phase(4,i);
    model4_chsip(find(model4(:,3)==60),4)=phase(3,i);
    model4_chsip(find(model4(:,3)==60),3)=amplitude(3,i);
    model4_chsip(find(model4(:,3)==500),4)=phase(2,i);
    model4_chsip(find(model4(:,3)==500),3)=amplitude(2,i);
    model4_chsip(find(model4(:,3)==90),4)=phase(1,i);
    model4_chsip(find(model4(:,3)==90),3)=amplitude(1,i);
    save(['model4_chsip_',num2str(f(i)),'.mod'],'model4_chsip','-ascii');
end
%% 


for i=1:length(f)
    forward=load(['forward_',num2str(f(i)),'.d']);
    c= [ 59 288 295 296 301 302 306 319 320 328 329 336 359 363 367 368 373 378 379 385 392 396 401 402 406 407 411 464 466 467 468 469 470 471 472 473 474 475 478 479 483];  
    ind=1:length(forward(:,1));
    c=ind(~ismember(ind,c));
    forward=forward(c,:);
    forward1=forward;
    neg1=find(forward(:,9)<0);
    neg2=find(forward(:,10)<0);

    %
    n=length(forward(:,1));
    amptmp=(forward(:,9).^2+forward(:,10).^2).^0.5;
    ampn=randn(n,1)*median(amptmp.*geofac)*0.05;

    phitmp=atan(forward(:,10)./forward(:,9))*1000;
    phin=randn(n,1);
    tmptmp=((amptmp.*geofac+ampn)./geofac).*exp((phitmp+phin)/1000*1i);
    forward(:,9)=real(tmptmp);
    forward(:,10)=imag(tmptmp);

    b=find(forward(:,10)<0);
    m2=length(b);
    while (m2~=0)
       amptmp=(forward1(b,9).^2+forward1(b,10).^2).^0.5;
       ampn=randn(m2,1)*median(amptmp.*geofac(b))*0.05;

       phitmp=atan(forward1(b,10)./forward1(b,9))*1000;
       phin=randn(m2,1);
       tmptmp=((amptmp.*geofac(b)+ampn)./geofac(b)).*exp((phitmp+phin)/1000*1i);
       forward(b,9)=real(tmptmp);
       forward(b,10)=imag(tmptmp);

       b=find(forward(:,10)<0);
       m2=length(b);
    end

    a=find(forward(:,9)<0);
    m=length(a)

%     forward(a,9)=(forward(:,9).*geofac+randn(n,1)*median(abs(forward(:,9).*geofac))*0.05)./geofac;
% 
%     a=find(forward(:,9)<0);
%     m=length(a);
%     while(m~=0)
%     forward(a,9)=(forward(a,9).*geofac(a)+randn(m,1)*median(abs(forward(:,9).*geofac))*0.05)./geofac(a);
%     a=find(forward(:,9)<0);
%     m=length(a);
%     end

%     forward(:,10)=(forward(:,10).*geofac+randn(n,1)*median(abs(forward(:,10).*geofac))*0.05)./geofac;
%     a=find(forward(:,10)<0);
%     m=length(a);
%     while(m~=0)
%     forward(a,10)=(forward(a,10).*geofac(a)+randn(m,1)*median(abs(forward(:,10).*geofac))*0.05)./geofac(a);
%     a=find(forward(:,10)<0);
%     m=length(a);
%     end

    save(['forward_chsipn_',num2str(f(i)),'_newnew.d'],'forward','-ascii');
    clear forward
end   %end loop over freq


%%
for i =1:length(f)
    forward=load(['forward_chsipn_',num2str(f(i)),'_newnew.d']);
    forward=forward([1:370 372:end],:);
    save(['forward_chsipn_',num2str(f(i)),'_newnew.d'],'forward','-ascii');
    clear forward
end
%%

R=0.20;

for i=1:length(f)
%   load(['inv_results_chsipn_',num2str(f(i)),'_newnew_1.5.mat']);
%   final_model(:,i)=conj(final.res_param1(:,final.itr))./(abs(final.res_param1(:,final.itr)).^2);

    final_model(:,i)=conj(final.d4_res_param1(:,i,final.itr))./(abs(final.d4_res_param1(:,i,final.itr)).^2);
    F(:,i)=0.01./(real(final_model(:,i))+5/R.*imag(final_model(:,i)));
end

% avgF=0.01./(real(final_model(:,length(f)))+5/R.*min(imag(final_model),[],2));
% 
avgF=mean(F(:,[1:4])')';
% avgF(find(avgF<0))=1;

%%
x=logspace(-3,2,100);
M=length(final_model(:,1));
D=3.8*1E-12;

%%
model=inline('-0.5*p(1).*cos(pi/2*(1-p(2)))./(cosh(p(2).*log(x*p(3)))+sin(pi/2*(1-p(2))))','p','x');
for i=1:M
    Y=imag(final_model(i,:));
    [beta3(i,:),r,J] =nlinfit2(w',Y',model,[-5/max(Y),0.5,1]);   
end
%%

%%
% ','Mn','c','tau0','w');
for i=1:M
    max_sigma_i(i)=feval(model,real(beta(i,:)),1/real(beta(i,3)));
end

max_sigma_i=feval(model,beta3,1./beta3(:,3));
avgF=0.01./(real(final_model(:,length(f)))+5/R.*max_sigma_i);
avgF(find(avgF<0))=1;
F(:,i)=0.01./(syn_sigma_inf+5/R.*imag(final_model(:,i)));
tau_0=(beta3(:,3));
tau_0(find(tau_0<0))=1e-3;
syn_k=D*tau_0/4./avgF;

%%
for i=1:M
    if(length(find(-imag(final_model(i,:))<0))>0)
        continue;
    end
    p=polyfit(log10(f),log10(-imag(final_model(i,:))),10);

%  p1=10.^polyval(p,log10(x));
%   loglog(f,-imag(final_model(i,:)),'o')
% hold on
% loglog(x,p1,'r')
% hold on
% title(num2str(i))

    k = polyder(p);
    kk = polyder(k);
    z=roots(k);
    t=z(find((polyval(kk,z)<0)&(imag(z)==0)));

%   if(max(t)<0.8) t=max(t); 

%   else
    t=t(find(t<2&t>-3));
%   [a,b]=max(polyval(p,t));
    t=min(t);%(b);
%   [a,b]=max(polyval(p,t));
%   if (t(b)<2 & t(b)>-3)
%      t=t(b);
%   else 
%      t=median(t(find(t>-3 & t<2)));
%   end
%   end

    if length(t)~=0;
       tau_0(i)=1/(2*pi*10^(t));
       syn_k(i)=D*tau_0(i)/4/avgF(i);
    else
       val_i=i
    end
end

%%
p1=10.^polyval(p,log10(x));
loglog(f,-imag(final_model(i,:)),'o')
hold on
loglog(x,p1,'r')

%%
title('one cell in formation 4')
xlabel('frequency (Hz)')
ylabel('quadruature conductivity (S/m)')

%%
real_k=[4.75*1e-13,2.11*1e-14,1.27*1e-15,5.9*1e-12];
model4_chsip(find(model4_chsip(:,3)==1000),4)=real_k(4);
model4_chsip(find(model4_chsip(:,3)==60),4)=real_k(3);
model4_chsip(find(model4_chsip(:,3)==500),4)=real_k(2);
model4_chsip(find(model4_chsip(:,3)==90),4)=real_k(1);

imagesc(final.param_x,final.param_y',log10(reshape(model4_chsip(:,4)',15,25)'))

