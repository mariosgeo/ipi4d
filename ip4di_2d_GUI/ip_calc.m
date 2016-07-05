function [input,mesh]=ip_calc(input,mesh)
disp ('------IP IVERSION STARTS-------');



% Creat new data set
input.real_data=input.real_data./(1-input.ip_data);


%Keep resistivity model
mesh.res_final_param=mesh.res_param2;

% Find new mean_res and update res_param1 and res_param2
mesh.mean_res=10^ (  mean (log10  (input.real_data(:,1)  )) );
mesh.res_param1(1:mesh.num_param)=mesh.mean_res;
mesh.res_param2(1:mesh.num_param)=mesh.mean_res;

% Update element properties
for i=1:mesh.num_param

    for j=1:mesh.num_elements

        if(mesh.icon(4,j)==abs(i))  ;mesh.prop(j)=mesh.res_param1(i); end
    end
end

end