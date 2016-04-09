%
% SUBROUTINE THAT DEFINES THE VARIABLE mesh.mean_res
%
%FL /!\ The value of mesh.mean_res must match homogeneous background, otherwise the solution is wrong.
% In previous version of Marios' code, mesh.mean_res was defined with respect to the data (apparent resisitivities) read in acqui file.
% Now we define it as the actual background value of the input model.

function mesh=define_mean_res(input,mesh)

 % read model from input file
 model=load(input.file_model_in);

 if input.cmplx_format==1
 %  input file contains real and imaginary part of resistivity
    rho_r=model(:,3);
    rho_i=model(:,4);
    mesh.model_in = complex (rho_r , rho_i);

 elseif input.cmplx_format==2
 % else, input file contains amplitude and phase of resistivity (in mrad)
    mag=model(:,3);
    phi=model(:,4);
    mesh.model_in = mag.*complex( cos(phi/1000) , sin(phi/1000) );
 end

 %define mean value
 if input.mean_res_flag==1
    %consider background value = value in 1st model cell
    mesh.mean_res=mesh.model_in(1);
 else
    %consider mean value
    mesh.mean_res=(mean(10.^(log10(mesh.model_in))));
 end


 % output mean value
 val_mean_res=mesh.mean_res


 %FL: lines below were initially in create_mesh3.m
 %FL /!\ Default values must match homogeneous background,
 %       otherwise solution is wrong
 mesh.init_res_param(1:mesh.num_param)=mesh.mean_res;
 mesh.res_param1(1:mesh.num_param,1)=mesh.mean_res;
 mesh.res_param2(1:mesh.num_param,1)=mesh.mean_res;
 mesh.prop(1:mesh.num_elements)=mesh.mean_res;

 if input.time_lapse_flag==1
     mesh.d4_res_param1=mesh.mean_res*ones(mesh.num_param,input.num_files);
     mesh.d4_res_param2=mesh.mean_res*ones(mesh.num_param,input.num_files);   
 end

end   %end function define_mean_res

