function [fem,mesh]=d4_prepare_data(itr,input,mesh,fem)





    fem.A=zeros(input.num_files*input.num_mes,input.num_files*mesh.num_param);
    fem.CC=zeros(input.num_files*mesh.num_param,input.num_files*mesh.num_param);
    fem.e=zeros(input.num_files*input.num_mes,1);
    fem.d4_array_model_data=zeros(input.num_mes,input.num_files);



h=waitbar(0,'4D calculations...');
    for j=1:input.num_files       
              
       waitbar(j/input.num_files,h);
        input.real_data=input.d4_real_data(:,j);

        % Update the prop matrix with each
%         for i=1:mesh.num_param
%             for k=1:mesh.num_elements
%                 if(mesh.elem_param(k)==abs(i))  ;mesh.prop(k)=mesh.d4_res_param1(i,j); end
%             end
%         end
        
        for i=1:mesh.num_param
            ind= mesh.elem_param==i;
            mesh.prop(ind)=mesh.d4_res_param1(i,j);
        end
        
        mesh.res_param1=mesh.d4_res_param1(:,j);
        [fem,ex]=mes_control_fast(itr,input,mesh,fem,j);
        
        % Keep each array_model_data and array_jacobian
        fem.d4_array_model_data(:,j)=fem.array_model_data;
          
        fem.e( (j-1)*input.num_mes+1 : j*input.num_mes,1)=log10(input.d4_real_data(:,j))-log10(fem.array_model_data);     
        fem.A( (j-1)*input.num_mes+1 : j*input.num_mes, (j-1)*mesh.num_param+1:j*mesh.num_param)=fem.array_jacobian;     
        if input.acb_flag==1
            [ctc,L1]=acb3d(input,mesh,fem);
        else
            ctc=input.lagrn*eye(mesh.num_param);            
        end
        fem.CC( (j-1)*mesh.num_param+1 : j*mesh.num_param, (j-1)*mesh.num_param+1:j*mesh.num_param)= ctc;

     end
        
close(h)

end