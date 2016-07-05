%Simple complex gradient code for complex valued A,b
function [x] = cgcomplex(A,b)
%initialize variables
L=length(A);
x = zeros(L,1); %x is your model vector; initial guess of zeros

%Jacobi preconditioner
P=speye(L);
DA=diag(A); %the actual preconditioner; in solution, solving (DA)^(-1)*r
T=sparse(L,L);
for ii=1:L
    if DA(ii)==0
        T(ii,ii)=0;
    else
        T(ii,ii)=(1/DA(ii))*P(ii,ii);   %inverse of preconditioner
    end
end

%set residuals
r = b-A*x;
z = T*r;
norm_b = sqrt(sum(b.*conj(b)));    %calculate the norm of the data
norm_r = sqrt(sum(r.*conj(r)));   %calculate the norm of the residual
err = norm_r/norm_b;      %divide to get the error between the two

rho0 = 1.0;             %set up rho0, used in finding the first beta
kk = 1;                  %counter

while (err >= 10e-23) && (kk < L)   %should find answer within N iterations
    
    rho = r.' * z;                  %calculate current rho
    if kk == 1
        beta = rho / rho0;          %for first iteration, initial search direction 
        p = z;                      %is along the residual, beta = rho
    else
        beta = rho / rhop;          %search direction updated using residual, beta, 
        p = z + beta * p;           %and previous search direction
    end

    q = A * p;
    alpha = rho / (p.' * q);         %set up minimizer of r'*A^-1*r
    x = x + (alpha * p);            %update model based on previous model, search direction
    r = r - (alpha * q);            %update residual based on previous residual
    z = T * r;
    rhop = rho;                     %reset rhop to use for new beta
    
    norm_r = sqrt(sum(r.*conj(r)));       %recalculate norm of the residual
    err = norm_r/norm_b;       %recalculate error
    
    kk = kk + 1;                      %increase counter
    
end
if kk >=L
    fprintf('N data reached with no convergence. Better luck next time.\n');
else
    fprintf('Convergence reached in %i iterations.\n',kk);
end

end