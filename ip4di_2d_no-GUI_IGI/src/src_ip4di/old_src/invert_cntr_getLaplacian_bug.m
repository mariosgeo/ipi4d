%
% FUNCTION INVERT_CNTR: UPDATE THE MODEL PARAMETERS FOLLOWING THE SELECTED OPTIMIZATION SCHEME
% Author: Marios Karaoulis, Colorado School of Mines

% Modified by Francois Lavoue', Colorado School of Mines, September 2015.
% Modifications include:
% 1) some explicit comments for more readability,
% 2) image-guided parameter update using the directional Laplacian filters defined in initcm2.m,
% 3) some change of notations, in particular:
%    - 'dx1' is now denoted 'Hx1' since it corresponds to the Hessian of the (total) objective function,
%    - 'dm1' is now denoted 'Gm1' because it corresponds to the gradient of the objective function,
%    - 'tmpmtx'  is now denoted 'Gm_mod'  because it is the model term of the gradient
%    - 'tmpmtx2' is now denoted 'Jm_term' because it is the J'm term used in Occam inversion

function [mesh,fem,input]=invert_cntr(itr,ip_cnt,input,mesh,fem)

%    /* ******************************************************** */
%    /* *************    Error weighting	******************** */
%    /* ******************************************************** */
% 
%    /* perform error weighting only when inv_flags are 3 or 4 */
%
% wd=abs((log10(input.real_data))-log10((fem.array_model_data)));
% % 
% % for i=1:input.num_mes
% %     if 100*abs(input.real_data(i)-fem.array_model_data(i))/input.real_data(i) <fem.rms_sum1
% %         wd(i)=fem.rms_sum1/100;
% %     end
% % end
% % 
% input.Wd=diag(1./wd);


 if input.inv_flag==3 || input.inv_flag==4
 %if Levenberg-Marquardt or Occam-Difference inversion: redefine data weighting matrix Wd

    sum1=0;
    sum2=0;

    %/* find initially sums */
    for i=1:input.num_mes
        % S1 = weighted sum of log residuals
        sum1 = sum1 + (input.Wd(i,i).^0.5)   *abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
        % S2 = weighted sum of sqrt( log residuals )
        sum2 = sum2 + (input.Wd(i,i).^0.25) *(abs(log10(input.real_data(i))-log10(fem.array_model_data(i))).^0.5);
    end

    hit=0;
    w_trial=zeros(input.num_mes);

    % define new (trial) data weighting matrix w_trial
    for i=1:input.num_mes
        w_trial(i,i) = (input.Wd(i,i)^0.5) / abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
        w_trial(i,i) = w_trial(i,i)*(sum1/sum2);   % w_trial ~ Wd^0.75 / abs( log(dobs)-log(dcal) )^0.5
        if w_trial(i,i)>input.Wd(i,i)    % keep min( w_trial , Wd ) (ie max. uncertainty)
           w_trial(i,i)=input.Wd(i,i);
           hit=hit+1;
        end
    end

    % /* find L1 ratio */
    sum1=0; sum2=0;
    for i=1:input.num_mes
        sum1 = sum1 + input.Wd(i,i)*abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
        sum2 = sum2 +  w_trial(i,i)*abs(log10(input.real_data(i))-log10(fem.array_model_data(i)));
    end

    % l1 = ratio between weighted sums of log residuals
    l1=sum1/sum2;

    %/* accept changes only when L1 ratio >1 */
    %   ie if the sum of log residuals weighted with input.Wd is greater than the one using w_trial
    %   ie keep the weights that provides the minimal weighted sum of log residuals...
    if (l1>1)
       input.Wd=w_trial;
    end
        
 end   %end if Levenberg-Marquardt or Occam-Difference inversion


 % 1st it.: init. Lagrangian multiplier
 if (ip_cnt==1 && itr==1); input.original_lagrn=input.lagrn; end
 if (ip_cnt==2 && itr==1); input.lagrn=input.original_lagrn; end


 % other it.: reduce Lagrangian multiplier (ie prior weight)
 if (itr~=1 && itr<5 && input.decrease_lagrn_flag==1) input.lagrn=input.lagrn/input.lagrn_reduction; end


 % Quasi Newton: update Jacobian
 if (itr>=2 && input.jacobian_flag==1)
    fem=quasi_newton(input,mesh,fem);
 end


 % compute Gauss-Newton Hessian (data part only)
 JTJ1=fem.array_jacobian.'*input.Wd*fem.array_jacobian;


 %------------------------------------------------------------------------------%
 % %/* find jacobian scale */
 %  sum=0;
 %  	for j=1:mesh.num_param
 % 	  for i=1:mesh.num_param
 % 		  tmp11=JTJ1(i,j);
 % 		  sum=sum+tmp11*tmp11;
 %       end
 %     end   
 %       
 %  jacscale1=sqrt(sum)
 % 
 %------------------------------------------------------------------------------%


 %------------------------------------------------------------------------------%
 %                                                                              %
 %                    COMPUTE GAUSS-NEWTON HESSIAN                              %
 %                                                                              %
 %------------------------------------------------------------------------------%
 %                                                                              %
 % Compute Hessian Hx1 = (J'J + lambda*Cm^-1) in                                %
 %                                                                              %
 % dm = - ( J'J + lambda*Cm^-1 )^-1 * (-J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m )
 %    = - (        Hx1         )^-1 * (                Gm1                      )
 %    = - (      Hessian       )^-1 * (              Gradient                   )
 %                                                                              %
 % (see Tarantola, 2005, eq. 3.51 and Karaoulis et al., eq. 19)                 %
 %                                                                              %
 % Reminder: ctc = Cm^-1                                                        %
 %                                                                              %
 %------------------------------------------------------------------------------%

 if input.acb_flag==0
 %  tmp_lagr=jacscale*input.lagrn;
    tmp_lagr=input.lagrn;
    Hx1=(JTJ1 + tmp_lagr*mesh.ctc); 
    ctc=mesh.ctc;

 elseif (input.acb_flag==1) 
 % if ACB is used, we recompute the covariance ctc,
 % locally weighted by the Lagrangian distribution L1
    [ctc,L1]=acb(input,mesh,fem);
    Hx1=(JTJ1+ctc); 

    % keep ACB's Lagrangian distribution L1 (vector of size mesh.num_param)
    % for further use with directional Laplacian filters
    fem.L1=L1;
 end
 %------------------------------------------------------------------------------%


 %------------------------------------------------------------------------------%
 %                                                                              %
 %                  COMPUTE MODEL TERM OF GRADIENT                              %
 %                                                                              %
 %------------------------------------------------------------------------------%
 %                                                                              %
 % Compute model term of Gradient Gm_mod = (lambda*Cm^-1*m) in                  %
 %                                                                              %
 % dm = - ( J'J + lambda*Cm^-1 )^-1 * (-J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m )
 %    = - (        Hx1         )^-1 * (                Gm1                      )
 %    = - (      Hessian       )^-1 * (              Gradient                   )
 %                                                                              %
 % (see Tarantola, 2005, eq. 3.51 and Karaoulis et al., eq. 19)                 %
 %                                                                              %
 % Reminder: ctc = Cm^-1                                                        %
 %                                                                              %
 %------------------------------------------------------------------------------%

 if input.inv_flag==2 
 % Gauss-Newton scheme
    if input.acb_flag==0
       Gm_mod=tmp_lagr*mesh.ctc*log10(mesh.res_param1);
    else
       Gm_mod=ctc*log10(mesh.res_param1);
    end
    Jm_term=zeros(input.num_mes,1);
 end


 if input.inv_flag==0 || input.inv_flag==5 
 % ?? or GN-Difference scheme
    Gm_mod=zeros(mesh.num_param,1);
    Jm_term=zeros(input.num_mes,1);
 end


 % Difference inversion
 if input.inv_flag==6 
    if input.acb_flag==0
       if itr>1 
          Gm_mod=tmp_lagr*ctc*(log10(mesh.res_param1) - log10(mesh.bgr_param));
       else
          Gm_mod=tmp_lagr*ctc*(log10(mesh.res_param1));
       end
    else   %acb_flag==1
       if itr>1
          Gm_mod=ctc*(log10(mesh.res_param1) - log10(mesh.bgr_param));
       else
          Gm_mod=ctc*(log10(mesh.res_param1));
       end
    end     

    Jm_term=zeros(input.num_mes,1);
 end


 if input.inv_flag==1
 % Occam scheme: model term is only J'm
 % (see Karaoulis et al., 2013, eq. 17)
    Gm_mod=zeros(mesh.num_param,1);
    Jm_term=fem.array_jacobian*log10(mesh.res_param1);
 end
 %------------------------------------------------------------------------------%


 %------------------------------------------------------------------------------%
 %                                                                              %
 %                         COMPUTE GRADIENT                                     %
 %                                                                              %
 %------------------------------------------------------------------------------%
 %                                                                              %
 % Compute Gradient Gm1 = (-J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m ) in        %
 %                                                                              %
 % dm = - ( J'J + lambda*Cm^-1 )^-1 * (-J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m )
 %    = - (        Hx1         )^-1 * (                Gm1                      )
 %    = - (      Hessian       )^-1 * (              Gradient                   )
 %                                                                              %
 % (see Tarantola, 2005, eq. 3.51 and Karaoulis et al., eq. 19)                 %
 %                                                                              %
 % Reminder: ctc = Cm^-1                                                        %
 %                                                                              %
 %------------------------------------------------------------------------------%

 if input.inv_flag~=5 && input.inv_flag~=6
 % if NON Difference inversion

    %     if itr>0
    %        a=imag(fem.array_model_data);
    %        b=find(a<0);
    %        c=find(a(b-1)<0);
    %        while(~isempty(c))
    %            b(c)=b(c)-1;
    %            c=find(a(b-1)<0);
    %        end
    %        
    %        a(find(a<0))=log10(a(b-1));
    %        Gm1=fem.array_jacobian.'*input.Wd*(complex(log10(real(input.real_data))-log10(real(fem.array_model_data)),...
    %         log10(imag(input.real_data)) -a) +Jm_term ) - Gm_mod;
    % % 10^(log10(abs(mesh.mean_res)))*exp(10^(log10(angle(mesh.mean_res))));
    %     else
    %             Gm1=fem.array_jacobian.'*input.Wd*(complex(log10(real(input.real_data))-log10(real(fem.array_model_data)),...
    %         log10(imag(input.real_data)) ) +Jm_term ) - Gm_mod;
    %    end

    %     !BE CAREFUL TO MINUS SIGN                                                      Occam-only     non-Occam
    % G = -    J'              *  Cd^-1 *(      d_obs           -      d_cal                 +   J'm  ) + lambda*Cm^-1*m
    Gm1 = -fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data) +Jm_term ) + Gm_mod;

    %Gm1=fem.array_jacobian.'*input.Wd*(log10(real(input.real_data))-log10(fem.array_model_data) +Jm_term ) - Gm_mod;   %FL: why "-G_mod"???
    %     a=imag(fem.array_model_data);
    %     b=find(a<0);
    %     c=find(a(b-1)<0);
    %     while(~isempty(c))
    %         b(c)=b(c)-1;
    %         c=find(a(b-1)<0);
    %     end
    %     
    %     a(find(a<0))=log10(a(b-1));
    %
    %     dm2=fem2.array_jacobian.'*input.Wd*(log10(imag(input.real_data))-log10(fem2.array_model_data)+Jm_term ) - Gm_mod;

 else
 % Difference inversion
     if itr>1
        % G = -    J'              *  Cd^-1 *(     (    d_obs       -      d_obs0          )-(     d_cal                -      d_cal0                    )
        Gm1 = -fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(input.real_data0)-log10(fem.array_model_data)+log10(input.array_model_data_bgr)) ...
            + Gm_mod;   %FL: "-G_mod" changed into "+G_mod"
        %   + lambda*Cm^-1*m
     else
        % G = -          J'        *  Cd^-1 *(      d_obs           -      d_cal                ) + lambda*Cm^-1*m
        Gm1 = -fem.array_jacobian.'*input.Wd*(log10(input.real_data)-log10(fem.array_model_data)) + Gm_mod;   %FL: "-G_mod" changed into "+G_mod"
     end
 end   %end if Difference inversion


 %------------------------------------------------------------------------------%
 %                                                                              %
 %                        COMPUTE MODEL UPDATE                                  %
 %                                                                              %
 %------------------------------------------------------------------------------%
 %                                                                              %
 % dm = - ( J'J + lambda*Cm^-1 )^-1 * (-J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m )
 %    = - (        Hx1         )^-1 * (                Gm1                      )
 %    = - (      Hessian       )^-1 * (              Gradient                   )
 %                                                                              %
 % (see Tarantola, 2005, eq. 3.51 and Karaoulis et al., eq. 19)                 %
 %                                                                              %
 %------------------------------------------------------------------------------%

 if input.image_guidance_flag~=2

 dx1 = -Hx1\Gm1;   % Solve H*dm = -G
 fem.dx1 = dx1;   % Keep this in case of Quasi Newton update

 %dx2=dx2\dm2;
 %fem.dx2=dx2;

 %%FL: comment this because now we compute explicitly the inverse covariance matrix mesh.ctc=Cm^-1=G'DG,
 %     so we can compute directly the model update
 %if(0)
 elseif input.image_guidance_flag==2
 % USE IMAGE-GUIDANCE THROUGH LOCAL LAPLACIAN
    %JZ: How to use this is after initializing (calling initcm2 before inversion)
    %    instead of doing in invert_cntr.m
    %       dX = (J'J + beta*Cm)^(-1) * J'*Cd*[G(X)-d)]
    %    do this:
    %       dX = reshape(interp.Sm.getLaplacianI(mesh.et2,reshape(J'*Cd*[G(X)-d)],mesh.m1,mesh.m2)',(J'J)',mesh.scale,mesh.s,beta)',mesh.num_param,1);
    % Here mesh.scale is how you want to scale the tensors, the larger the more emphasize on structure.
    %    FL: this is redundant with the role of the Lagrangian multiplier beta=input.lagrn, better don't use mesh.scale.
    %    FL2JZ: can you confirm my understanding?
    % mesh.s is using structure-oriented semblance s2 to (locally) scale the tensors, if don't wanna use can just put [] instead of mesh.s
    %    FL: this is redundant with the spatially-varying Lagrangian defined by ACB (but slightly different, though...)
    %    FL2JZ: can you confirm my understanding?

    %define scale and beta parameters
    %such as to avoid redundancy...
    scale=1;
    beta=input.lagrn;

    Gm1 = reshape(Gm1,mesh.m1,mesh.m2);   %reshape gradient in matrix form
    %dm = ( J'J + lambda*Cm^-1 )^-1 * ( J'*Cd^-1*[d_obs-d_cal] + lambda*Cm^-1*m )
    %   = ( J'J + lambda*Cm^-1 )^-1 * (                Gm1                      )
    %   = (      Hessian       )^-1 * (              Gradient                   )
    dx1 = interp.Sm.getLaplacianI(mesh.et2,Gm1',JTJ1',scale,[],beta);
    dx1 = reshape(dx1',mesh.num_param,1);
 %end   %end if(0) to comment block of lines

 end   %end if image guidance
 %------------------------------------------------------------------------------%


 %------------------------------------------------------------------------------%
 %                                                                              %
 %                       UPDATE MODEL PARAMETERS
 %                                                                              %
 %------------------------------------------------------------------------------%

 for i=1:mesh.num_param
     if input.inv_flag==2 ||input.inv_flag==6 || input.inv_flag==0 || input.inv_flag==5
     %if GN, GN-Difference, or ??
        b=1; 
        a=10^(log10(mesh.res_param2(i)) + b*dx1(i));

        if imag(a)>0
        %if imag. part of updated parameter is >0, then simply update
           mesh.res_param1(i)=a;
        else 
        %if imag. part of updated parameter is <0,
        %then find a descent step b such that imag(m+b*dm)>0
           while (imag(a)<0)
              b=b/2; 
              a=10^(log10(mesh.res_param2(i)) + b*dx1(i));
           end
           mesh.res_param1(i)=a;
        end

% % if itr>0
% %  mesh.res_param1(i)=complex(10^(log10(real(mesh.res_param2(i)))+real(dx1(i))),10^(log10(imag(mesh.res_param2(i)))+imag(dx1(i))));
% % else
% %     mesh.res_param1(i)=complex(10^(log10(abs(mesh.res_param2(i)))+real(dx1(i))),10^(imag(dx1(i))));
% % end
%     mesh.res_param1(i)=10^(log10(mesh.res_param2(i)) + dx1(i));
%   mesh.res_param1(i)=complex(10^(log10(real(mesh.res_param2(i))) + dx1(i)),10^(log10(imag(mesh.res_param2(i))) + dx2(i)));
%       mesh.res_param1(i)=10^log10(complex(10.^(log10(real(input.real_data))),...
%         10.^(log10(imag(input.real_data))))

     else
     % if other optimization than GN, GN-Diff or ??
        mesh.res_param1(i)=10^(dx1(i));
     end   %end if optimization type

%     if imag(mesh.res_param1(i))>0 ; mesh.res_param1(i)=complex(real(mesh.res_param1(i)),-0.01); end

     if input.limit_res==1
     % apply rough bound constraints
        if mesh.res_param1(i)>input.max_res mesh.res_param1(i)=input.max_res; end
        if mesh.res_param1(i)<input.min_res mesh.res_param1(i)=input.min_res; end
     end

 end   %end for num_param


 % compute resolution matrix
 % R = ( Hd + lambda*Hm )^-1 * Hd
 %   = (J'J + lbda*Cm^-1)^-1 * J'J
 fem.resolution_matrix = Hx1\JTJ1;   % solve H*R=J'J, ie R=(H^-1)*J'J

 % data resolution matrix (without model term, just add eps*Id to stabilize inversion)
 % R = ( Hd + eps*I )^-1 * Hd
 %   = (J'J + eps*I )^-1 * J'J
 fem.resolution_matrix_data = (JTJ1 + 0.005*eye(mesh.num_param))\JTJ1;


 sprintf('**  ITERATION =>  %d  **\n',itr)
 %tt=sprintf('**  ITERATION =>  %d  **\n',itr);
 %tt=cellstr(tt);
 %drawnow;          

end   %end function

