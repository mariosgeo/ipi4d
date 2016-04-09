%
% FUNCTION KIM_INVERSION2: update model parameters in case of time-lapse inversion
%                          (equivalent of invert_cntr.m in case of non-time-lapse).
%
% Author: Marios Karaoulis, Colorado School of Mines

% Modified by Francois Lavoue', Colorado School of Mines, October 2015.
% Modifications include:
% 1) some comments for more readability,
% 2) the consideration of the model term Cm^-1*log(m) in the gradient.
% 3) the introduction of the spatially-varying Lagrangian weight through ACB (TO DO).

function [mesh]=kim_inversion2(input,mesh,fem)

% Lagrangian multiplier associated to time-lapse regularisation
alpha=input.gamma;

% init. time-lapse models (current+updated)
U=zeros(input.num_files*mesh.num_param,1);
U2=zeros(input.num_files*mesh.num_param,1);


% recast current time-lapse models in a single vector
for i=1:input.num_files
    U( (i-1)*mesh.num_param+1 : i*mesh.num_param , 1 ) = mesh.d4_res_param1(:,i);
end


% compute model update
if input.atc_flag==0
   %MK's formulation (a model term has been forgotten in the gradient)
   %dm =-(     Hessien data      + Hm     + Hessien time-reg)^-1 * (     Gradient data     + Gradient time-regularisation )
   %DU = (fem.A.'*input.Wd*fem.A + fem.CC + alpha*mesh.M'*mesh.M)\ (fem.A.'*input.Wd*fem.e - alpha*mesh.M'*mesh.M*log10(U));

   %dm =-(     Hessien data      + Hm     + Hessien time-reg )^-1 * (    Gradient data       +   Gradient model           )
   %dm =-(     J'*      Wd*    J + Cm^-1  + gamma*      M'*M )^-1 * (   -J'*        Wd*  res +      lambda*Cm^-1*log(m) + gamma*     M'*M     *log(m)  )
    DU =-(fem.A.'*input.Wd*fem.A + fem.CC + alpha*mesh.M'*mesh.M) \ (-fem.A.'*input.Wd*fem.e + input.lagrn*fem.CC*log10(U) + alpha*mesh.M'*mesh.M*log10(U));

%FL DEBUG WARNING: SHOULD IT BE   gamma * M * log(m)   OR   gamma * M * log(m)   in the gradient???
%   +++ INTRODUIRE L'ACB (APPLIQUER fem.L1 A fem.CC)

elseif input.atc_flag==1
% use Active-Time Constraints

    % read ACT matrix
    ACT=load(input.par_nm);
    ACT=ACT.act;

   %dm =-(     Hessien data      + Hm     + Hessien ACT      )^-1 * (   Gradient data   +   Gradient ACT             )
   %dm = (     J'*      Wd*    J + Cm^-1  +      M'*ACT*M    )^-1 * (J' *      Wd*  res -      M'*ACT*M     *log(m)  )
    DU = (fem.A.'*input.Wd*fem.A + fem.CC + mesh.M'*ACT*mesh.M)\(fem.A.'*input.Wd*fem.e - mesh.M'*ACT*mesh.M*log10(U));
end


for i=1:input.num_files*mesh.num_param
%   U2(i,1)=10^(log10(U(i)) + DU(i));
    b=1; 
    a=10^(log10(U(i)) + b*DU(i));
    if imag(a)>0
       U2(i,1)=a;

    else 
       while (imag(a)<0)
          b=b/2; 
          a=10^(log10(U(i)) + b*DU(i));
       end
       U2(i,1)=a;
    end
end


for k=1:input.num_files
    mesh.d4_res_param1(:,k)= U2 ( (k-1)*mesh.num_param+1:k*mesh.num_param,1);   
end


% 
% for k=1:input.num_files
%     for i=1:mesh.num_param
%         if imag(mesh.d4_res_param1(i,k))>0 
%             mesh.d4_res_param1(i,k)=complex(real(mesh.d4_res_param1(i,k)),-0.01);
%         end
%     end
% end
% 

end   %end function

