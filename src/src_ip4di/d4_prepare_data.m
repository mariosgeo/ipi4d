%
%            COMPUTES SYNTHETIC DATA 
% FOR EACH DATA SET IN THE TIME-LAPSE INVERSION
%
function [fem,mesh]=d4_prepare_data(itr,input,mesh,fem)

    % init. arrays for Jacobians, model covariance, residuals and data, respectively
    fem.A=zeros(input.num_files*input.num_mes,input.num_files*mesh.num_param);
    fem.CC=zeros(input.num_files*mesh.num_param,input.num_files*mesh.num_param);
    fem.e=zeros(input.num_files*input.num_mes,1);
    fem.d4_array_model_data=zeros(input.num_mes,input.num_files);

    % loop over data sets
    for j=1:input.num_files       

        % observed data
        input.real_data=input.d4_real_data(:,j);

        % update mesh properties with updated model parameters
        for i=1:mesh.num_param
            ind= mesh.icon(4,:)==i;
            mesh.prop(ind)=mesh.d4_res_param1(i,j);
        end
        mesh.res_param1=mesh.d4_res_param1(:,j);

        % SOLVE FORWARD PROBLEM
        fem=mes_control_fast(itr,input,mesh,fem,j);

        % store synthetic data
        fem.d4_array_model_data(:,j)=fem.array_model_data;

        % store log residuals
        fem.e( (j-1)*input.num_mes+1 : j*input.num_mes, 1 ) = log10(input.d4_real_data(:,j))-log10(fem.array_model_data);     

        % store Jacobian matrices
        fem.A( (j-1)*input.num_mes+1 : j*input.num_mes, (j-1)*mesh.num_param+1 : j*mesh.num_param ) = fem.array_jacobian;     

        % define model covariance matrix
        if input.acb_flag==1
        % ACB: define a locally-varying Lagrangian distribution based on data resolution
        %      and apply it to model covariance matrix
           input.Wd = input.Wd_d4( (j-1)*input.num_mes+1 : j*input.num_mes , (j-1)*input.num_mes+1 : j*input.num_mes );
           [ctc,fem.L1]=acb(input,mesh,fem);
        else
        % keep the model covariance matrix as it is
        %NB: do NOT include the scalar Lagrange multiplier here,
        %    it will be applied later on, in kim_inversion2.m
            %ctc=input.lagrn*eye(mesh.num_param);
            ctc=mesh.ctc;
        end

        % store time-lapse model covariance matrix
        fem.CC( (j-1)*mesh.num_param+1 : j*mesh.num_param , (j-1)*mesh.num_param+1 : j*mesh.num_param ) = ctc;

     end

end   %end function
