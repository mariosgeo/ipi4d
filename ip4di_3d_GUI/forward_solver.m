%The main idea is to solve the forward proble two times. One for normalization and one for calculations 


function fem=forward_solver(input,mesh)

fem=[]; % initialize fem as empty matrix
final=[];
if matlabpool('size') == 0
    matlabpool open
end
input.jacobian_flag=0; 

mesh.prop(:)=mesh.mean_res;


[fem,ex]=mes_control_fast(1,input,mesh,fem,0);



%for i=1:mesh.num_param
%
%    for j=1:mesh.num_elements
%
%        if(mesh.elem_param(j)==abs(i))  ;mesh.prop(j)=mesh.bgr_param(i); end
%    end
%end
 


    for i=1:mesh.num_param
        ind= mesh.elem_param==i;
        mesh.prop(ind)=mesh.bgr_param(i);

    end
%solve again for the user defined resistivity
input.jacobian_flag=0;
[fem,ex]=mes_control_fast(2,input,mesh,fem,0);
   
if input.sip_flag==0    
final=[input.ax input.ay input.az input.bx input.by input.bz input.mx input.my input.mz input.nx input.ny input.nz fem.array_model_data];
else
final=[input.ax input.ay input.az input.bx input.by input.bz input.mx input.my input.mz input.nx input.ny input.nz real(fem.array_model_data) imag(fem.array_model_data)];
end

save('forward.txt','final','-ascii');
% final=[input.ax;input.ay;input.bx;input.by;input.mx;input.my;input.nx;input.ny;fem.array_model_data];
   
   
   
   

end