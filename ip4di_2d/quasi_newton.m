function fem=quasi_newton(input,mesh,fem)

    dd=0;
    for i=1:input.num_param
        dd=dd+fem.dx1(i)*fem.dx1(i);
    end

    dm=log10(fem.array_model_data) - log10(fem.oldmes) - fem.array_jacobian*fem.dx1;

    %array_jacobian=array_jacobian + dm*dx1/dd;
    for i=1:input.num_mes
        for j=1:mesh.num_param
            fem.array_jacobian(i,j)=fem.array_jacobian(i,j) + dm(i)*fem.dx1(j)/dd;
        end
    end
    fem.oldmes=fem.array_model_data;


end