function [mesh]=kim_inversion2(input,mesh,fem)
%clear A CC U e M DU U2


alpha=0.01;


U=zeros(input.num_files*mesh.num_param,1);
U2=zeros(input.num_files*mesh.num_param,1);

for i=1:input.num_files
    U((i-1)*mesh.num_param+1 : i*mesh.num_param,1)=mesh.d4_res_param1(:,i);
end


if input.atc_flag==0
    DU=(fem.A.'*input.Wd*fem.A + fem.CC + alpha* mesh.M'*mesh.M)\(fem.A.'*input.Wd*fem.e - alpha*mesh.M'*mesh.M*log10(U));
elseif input.atc_flag==1
    % Here read ACT matrix?
    ACT=load(input.par_nm);
    ACT=ACT.act;
    DU=(fem.A.'*input.Wd*fem.A + fem.CC + mesh.M'*ACT*mesh.M)\(fem.A.'*input.Wd*fem.e - mesh.M'*ACT*mesh.M*log10(U));
end


for i=1:input.num_files*mesh.num_param
    U2(i,1)=10^(log10(U(i)) + DU(i));
end


for k=1:input.num_files
  mesh.d4_res_param1(:,k)= U2 ( (k-1)*mesh.num_param+1:k*mesh.num_param,1);   
end








for k=1:input.num_files
    for i=1:mesh.num_param
        if imag(mesh.d4_res_param1(i,k))>0 
            mesh.d4_res_param1(i,k)=complex(real(mesh.d4_res_param1(i,k)),-0.01);
        end
    end
end










end
