function [k,g,obj,err]=xu_inversion(rdis,k)
rdis=rdis'; 
%set the maximum number of iterations for the optimization routine
itsmax = 25;

%Max number of radii to search over
Max_num =2000;

%Number of linesearch steps
lsnum = 10;
%Line Search parameters
%lower bound
ls_low_lim = 0.01;
%upper bound
ls_up_lim =1;

num=length(k);



%K values matrix   
Km = ones(size(rdis,1),1) *(k(:))';
rdis2=rdis*ones(1,num);
%Identity vector
I = ones(size(rdis,1),1);


%Calculate the A matrix
A = rdis2.*real(besselk(0,rdis2.*Km));

%%Estimate g for the given K values (THIS IS NOT G)
v = A*((A'*A)\(A'*I));
%Evaluate the objective function for the initial guess
obj(1) = (1-v)'*(1-v);



k0=k;
%Start counter and initialize the optimization
its = 1; %iteration counter
knew=k0; %updated k vector
stop =0; %Stopping toggle incase A becomes illconditioned
reduction = 1; %Variable for ensure sufficent decrease between iterations
                    % Optimization terminates if objective function is not
                    % reduced by at least 5% at each iteration

   while obj >1e-5 & its <itsmax & stop ==0 & reduction > 0.05; 
  
    %%Create the derivative matrix
    dvdk = zeros(length(v),num);
    for i = 1:num;
        Ktemp = Km;
        Ktemp(:,i) = 1.05*(Ktemp(:,i));
         %form an new A matrix
         A = rdis2.*real(besselk(0,rdis2.*Ktemp));
             L = A'*A;
        
        %%Estimate g for the given K values
       vT = A*((L)\(A'*I));
       
       %Calculate the derivative for the appropriate column
        dvdk(:,i) = (vT-v)./(Ktemp(:,i)-Km(:,i));
    end;
    
    %Apply some smallness regularization
    h = dvdk'*(I-v)+1e-8*eye(length(knew))*knew(:);
    dk = (dvdk'*dvdk+1e-8*eye(length(knew)))\h;
    
    %Perform a line-search to maximize the descent 
    for j =[1:lsnum];
        warning off;
        ls =linspace(ls_low_lim,ls_up_lim,lsnum);
    ktemp =knew(:) +ls(j)*dk(:);

    Km = ones(size(rdis,1),1) *(ktemp(:))';
    %Matrix of ones
    %Calculate the A matrix
    A = rdis2.*real(besselk(0,rdis2.*Km));
     L = A'*A;
    %%Estimate g for the given K values
  
    v = A*(L\(A'*I));
    objt = (1-v)'*(1-v);
    ls_res(j,:) = [objt ls(j)];
  
    warning on;
    end;
   
    %Find the smallest objective function from the line-search
    [b,c] = (min(ls_res(:,1)));
         
    %Create a new guess for k
    knew =knew(:) +ls(c)*dk(:);
    %eval obj funct
    Km = ones(size(rdis,1),1) *(knew(:))';
    
    %Calculate the A matrix
    A = rdis2.*real(besselk(0,rdis2.*Km));
    %%Estimate g for the given K values
      
    v = A*((A'*A)\(A'*I));
    obj(its+1)= (1-v)'*(1-v);
    reduction = obj(its)./obj(its+1)-1;         
    its = its+1;
    %Check the conditioning of the matrix
    if rcond(A'*A) < 1e-20;
        knew =knew(:) -ls(c)*dk(:);
        stop = 1;
    end;
end;

%Get the RMS fit 
err = sqrt(obj./length(rdis));
%The final k values
k =abs(knew);

Km = ones(size(rdis,1),1) *(k(:))';
%Reform A to obtian the final g values
A = rdis2.*real(besselk(0,rdis2.*Km));
%Calculate g values
g = ((A'*A)\(A'*I));
k=reshape(k,num,1);