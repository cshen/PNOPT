function Hx = LbfgsProd(sPrev, yPrev, et) 
% LbfgsProd : Product with L-BFGS Hessian approximation
% 
%   $Revision: 0.1.1 $  $Date: 2012/06/15 $
% 
  l = size(sPrev,2);
  L = zeros(l);
  for k = 1:l;
    L(k+1:l,k) = sPrev(:,k+1:l)'*yPrev(:,k);
  end
  d1 = sum(sPrev.*yPrev);
  d2 = sqrt(d1);
  
% %   Mark Schmidt's code (slow)
%   Qty1 = et*(sPrev'*sPrev);
%   Qty2 = [et*sPrev, yPrev];
%   Qty3 = [Qty1,  L;
%           L'  , -diag(d1) ];
%   Hx   = @(x) et*x - Qty2*(Qty3\(Qty2'*x));
%   
  R    = chol(et*(sPrev'*sPrev) + L*(diag(1./d1)*L'),'lower');
  R1   = [diag(d2)       zeros(l);
          -L*diag(1./d2) R];
  R2   = [-diag(d2)      diag(1./d2)*L';
          zeros(l)       R'];
  Qty2 = [yPrev, et*sPrev];
  Hx   = @(x) et*x - Qty2*(R2\(R1\(Qty2'*x)));
  