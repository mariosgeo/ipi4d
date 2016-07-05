%
% INIT. 4D MODEL
%
% IN PROGRESS........

function [input,mesh]=init_4d_model(input,mesh)

 % read model from input file
 % (same format as for 4d data)
 tmp=importdata(input.file_model_in);
 num_files=length(tmp);

 % check
 if num_files~=input.num_files
    disp(sprintf('NB OF INITIAL MODELS = %i =/= %i, ABORT.',num_files,input.num_files))
    return
 end

 % init. arrays
 mesh.mean_res_4d=zeros(input.num_files,1);
 mesh.model_in_4d=zeros(mesh.num_param,input.num_files);
 mesh.d4_res_param1=zeros(mesh.num_param,input.num_files);

 % recast models in a 3rd-order tensor [i_row=i_mod,i_col=real/imag,i_file]
 for i=1:input.num_files
    model=tmp{i};

    if input.cmplx_format==1
    %  input file contains real and imaginary part of resistivity
       rho_r=model(:,3);
       rho_i=model(:,4);
       mesh.model_in_4d(:,i) = complex (rho_r , rho_i);

    elseif input.cmplx_format==2
    % else, input file contains amplitude and phase of resistivity (in mrad)
       mag=model(:,3);
       phi=model(:,4);
       mesh.model_in_4d(:,i) = mag.*complex( cos(phi/1000) , sin(phi/1000) );
    end

    %define mean value
    if input.mean_res_flag==1
       %consider background value = value in 1st model cell
       mesh.mean_res_4d(i)=mesh.model_in_4d(1,i);
    else
       %consider mean value
       mesh.mean_res_4d(i)=mean(10.^(log10(mesh.model_in_4d(:,i))));
    end

    % init. model
    mesh.d4_res_param1(:,i)=mesh.mean_res_4d(i);

 end   %end for i=1:input.num_files

 %if input.debug==1
     % output mean value
     val_mean_res=mesh.mean_res_4d
 %end

 %FL: lines below were initially in create_mesh3.m
 %FL /!\ Default values must match homogeneous background,
 %       otherwise solution is wrong
 %mesh.res_param1(1:mesh.num_param,1)=mesh.mean_res;
 %mesh.res_param2(1:mesh.num_param,1)=mesh.mean_res;
 mesh.prop(1:mesh.num_elements)=mesh.mean_res_4d(1);
 mesh.d4_res_param2=mesh.d4_res_param1;

end   %end function define_mean_res

