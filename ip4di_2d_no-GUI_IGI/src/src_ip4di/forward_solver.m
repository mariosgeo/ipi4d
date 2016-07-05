%MK: The main idea is to solve the forward problem two times. One for normalization and one for calculations.

function fem=forward_solver(input,mesh)

fem=[];

% Calculate k and g
[mesh.k,mesh.g]=rwaven(mesh.probe_spacing,mesh.probe_spacing*mesh.max_n,4);

%assign constant value to all elements
mesh.prop(:)=mesh.mean_res;

%compute data in homogeneous medium
fem=mes_control_fast(0,input,mesh,fem,0);

%update now res_param
mesh.res_param1=input.bgr_res_param;
for i=1:mesh.num_param
    for j=1:mesh.num_elements
        if(mesh.icon(4,j)==abs(i)); mesh.prop(j)=mesh.res_param1(i); end
    end
end
 
%solve again for the user defined resistivity
fem=mes_control_fast(2,input,mesh,fem,0);

if input.sip_flag==0    
    final=[input.ax input.az input.bx input.bz input.mx input.mz input.nx input.nz real(fem.array_model_data)];
    model=[mesh.param_x mesh.param_y real(mesh.res_param1)];
else
    final=[input.ax input.az input.bx input.bz input.mx input.mz input.nx input.nz real(fem.array_model_data) imag(fem.array_model_data)];
    model=[mesh.param_x mesh.param_y real(mesh.res_param1)  imag(mesh.res_param1)];
end

save(input.file_data_out,'final','-ascii');
% final=[input.ax;input.ay;input.bx;input.by;input.mx;input.my;input.nx;input.ny;fem.array_model_data];
save(input.file_model_out,'model','-ascii');   

end   %end function forward_solver
