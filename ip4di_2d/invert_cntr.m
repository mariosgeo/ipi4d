function [mesh,fem,input]=invert_cntr(itr,tmp_flag,ip_cnt,input,mesh,fem)




%    /* ******************************************************** */
%    /* *************    Error weighting	******************** */
%    /* ******************************************************** */
% 
%      /* perform error weighting only when inv_flags are 3 or 4 */



% wd=abs((log10(input.real_data))-log10((fem.array_model_data)));
% % 
% % 
% % 
% % 
% % for i=1:input.num_mes
% %     if 100*abs(input.real_data(i)-fem.array_model_data(i))/input.real_data(i) <fem.rms_sum1
% %         wd(i)=fem.rms_sum1/100;
% %     end
% % end
% % 
% input.Wd=diag(1./wd);


if input.inv_flag==3 || input.inv_flag==4

sum1=0;
sum2=0;

%/* find initially sums */
          for i=1:input.num_mes
              sum1=sum1+ (input.Wd(i,i).^0.5)   *abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
              sum2=sum2+ (input.Wd(i,i).^0.25) *(abs(log10(input.real_data(i))-log10(fem.array_model_data(i))).^0.5);
          end
        hit=0;
        w_trial=zeros(input.num_mes);
        for i=1:input.num_mes
            w_trial(i,i)=(input.Wd(i,i)^0.5) / abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
            w_trial(i,i)=w_trial(i,i)*(sum1/sum2);
            if w_trial(i,i)>input.Wd(i,i) 
                w_trial(i,i)=input.Wd(i,i);
            hit=hit+1;
            end
        end

        % /* find L1 ratio */
             sum1=0; sum2=0;
             for i=1:input.num_mes
              sum1=sum1+input.Wd(i,i)*abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
              sum2=sum2+w_trial(i,i)*abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
             end
               l1=sum1/sum2;

        %/* accept changes only when L1 ratio >1 */
        if(l1>1)
        
            input.Wd=w_trial;
        
        end
        
end
        
        
        
        
        
        
        
        

 if (ip_cnt==1 && itr==1) ;input.original_lagrn=input.lagrn; end
 if (ip_cnt==2 && itr==1) ;input.lagrn=input.original_lagrn; end
 
 
 if (itr~=1 && itr<5 && tmp_flag==1) input.lagrn=input.lagrn/input.lagrn_reduction; end
 
 % Quasi Newton update Jacobian
 if(itr>=2 && input.jacobian_flag==1) ;fem=quasi_newton(input,mesh,fem); end

 
 
 JTJ=fem.array_jacobian.'*input.Wd*fem.array_jacobian;
 
 %				  /* find jacobian scale */
 sum=0;
 	for j=1:mesh.num_param
	  for i=1:mesh.num_param
		  tmp11=JTJ(i,j);
		  sum=sum+tmp11*tmp11;
      end
    end   
      
 jacscale=sqrt(sum);

 
 
  %--------------------ACB+dx1----------------------------------------------
if input.acb_flag==0
     tmp_lag=jacscale*input.lagrn;
     tmp_lag=input.lagrn;
     dx1=(JTJ + tmp_lag*mesh.ctc); 
     ctc=mesh.ctc;
elseif (input.acb_flag==1) 
    [ctc,L1]=acb(input,mesh,fem);
    dx1=(JTJ+ctc); 
    % keep ACB
    fem.L1=L1;
end
 


 %---------------C*dx------------------------------------------------------      
 if input.inv_flag==2 
     if input.acb_flag==0
        tmpmtx=tmp_lag*mesh.ctc*log10(mesh.res_param1);
     else
        tmpmtx=ctc*log10(mesh.res_param1);
     end
     tmpmtx2=zeros(input.num_mes,1);
 end
 
  if input.inv_flag==0 || input.inv_flag==5 
     tmpmtx=zeros(mesh.num_param,1);
     tmpmtx2=zeros(input.num_mes,1);
 end
 
 if input.inv_flag==6 
     
    
     
     if input.acb_flag==0
         if itr>1 
             tmpmtx=tmp_lag*ctc*(log10(mesh.res_param1) - log10(mesh.bgr_param));
         else
              tmpmtx=tmp_lag*ctc*(log10(mesh.res_param1));
         end
     else
        if itr>1
            tmpmtx=ctc*(log10(mesh.res_param1) - log10(mesh.bgr_param));
        else
            tmpmtx=ctc*(log10(mesh.res_param1));
        end
     end     
     
     tmpmtx2=zeros(input.num_mes,1);
 end
 
 
 if input.inv_flag==1 
     tmpmtx=zeros(mesh.num_param,1);
     tmpmtx2=fem.array_jacobian*log10(mesh.res_param1);
 end
%--------------------------------------------------------------------------

 
 
%---------------------Jt*dy------------------------------------------------
 if input.inv_flag~=5 && input.inv_flag~=6
     dm1=fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data) +tmpmtx2 ) - tmpmtx;
 else
     if itr>1 
         dm1=fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(input.real_data0)-log10(fem.array_model_data)+log10(input.array_model_data_bgr)) - tmpmtx;
     else
         dm1=fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data)) - tmpmtx;
     end
 end
 



dx1=dx1\dm1;
fem.dx1=dx1; % Keep this in case of Quasi Newton update




    for i=1:mesh.num_param
        if input.inv_flag==2 ||input.inv_flag==6 || input.inv_flag==0 ||input.inv_flag==5
            mesh.res_param1(i)=10^(log10(mesh.res_param2(i)) + dx1(i));
           
        else
              mesh.res_param1(i)=10^(dx1(i));
        end
%     if imag(mesh.res_param1(i))>0 ; mesh.res_param1(i)=complex(real(mesh.res_param1(i)),-0.01); end
        if input.limit_res==1
            if mesh.res_param1(i)>input.max_res mesh.res_param1(i)=input.max_res; end
            if mesh.res_param1(i)<input.min_res mesh.res_param1(i)=input.min_res; end
        end

    end


    
% keep resolution matrix
fem.resolution_matrix=dx1\(JTJ);



tt=sprintf('**  ITERATION =>  %d  **\n',itr);
tt=cellstr(tt);
drawnow;          



end











