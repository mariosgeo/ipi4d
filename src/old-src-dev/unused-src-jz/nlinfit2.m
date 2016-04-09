function [beta,r,J] = nlinfit2(X,y,model,beta0)
%NLINFIT Nonlinear least-squares data fitting by the Gauss-Newton method.
%   NLINFIT(X,Y,'MODEL',BETA0) finds the coefficients of the nonlinear 
%   function described in MODEL. MODEL is a user supplied function having 
%   the form y = f(beta,x). That is MODEL returns the predicted values of y
%   given initial parameter estimates, beta, and the independent variable, X.   
%   [BETA,R,J] = NLINFIT(X,Y,'MODEL',BETA0) returns the fitted coefficients
%   BETA the residuals, R, and the Jacobian, J, for use with NLINTOOL to
%   produce error estimates on predictions.

%   B.A. Jones 12-06-94.
%   Copyright (c) 1993-96 by The MathWorks, Inc.
%   $Revision: 2.7 $  $Date: 1996/09/29 17:33:27 $

n = length(y);
if min(size(y)) ~= 1
   error('Requires a vector second input argument.');
end
y = y(:);

p = length(beta0);
beta0 = beta0(:);

J = zeros(n,p);
beta = beta0;
betanew = beta + 1;
maxiter = 200;
iter = 0;
betatol = 1.0E-8;
rtol = 1.0E-8;
sse = 1;
sseold = sse;

while (norm((betanew-beta)./(beta+sqrt(eps))) > betatol | abs(sseold-sse)/(sse+sqrt(eps)) > rtol) & iter < maxiter
   if iter > 0, 
      beta = betanew;
   end

   iter = iter + 1;
   evalstr = model;
   yfit = feval(evalstr,beta,X);
   r = y - yfit;
   sseold = r'*r;

   for k = 1:p,
      delta = zeros(size(beta));
      delta(k) = sqrt(eps)*beta(k);
      evalstr1 = model;
      yplus =feval(evalstr1,beta+delta,X);
      J(:,k) = (yplus - yfit)/(sqrt(eps)*beta(k));
   end

   Jplus = [J;1*eye(p)]; %1.0E-2
   rplus = [r;zeros(p,1)];

   % Levenberg-Marquardt type adjustment 
   % Gauss-Newton step -> J\r
   % LM step -> inv(J'*J+constant*eye(p))*J'*r
   step = Jplus\rplus;
%    step=inv(Jplus'*Jplus+1000*eye(p))*Jplus'*rplus;
   betanew = beta + step;
   evalstrnew = model;
   yfitnew = feval(evalstrnew,betanew,X);
   rnew = y - yfitnew;
   sse = rnew'*rnew;
   iter1(iter) = 0;
   while sse > sseold & iter1(iter) < 12
      step = step/sqrt(10);
      betanew = beta + step;
      evalstrnew = model;
      yfitnew = feval(evalstrnew,betanew,X);
      rnew = y - yfitnew;
      sse = rnew'*rnew;
      iter1(iter) = iter1(iter) + 1;
   end
end
