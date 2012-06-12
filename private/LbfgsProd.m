function Hx = LbfgsProd(sPrev, yPrev, et) 
% LbfgsProd : Product with L-BFGS Hessian approximation
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $
% 
  l = size(sPrev,2);
  L = zeros(l);
  for k = 1:l;
    L(k+1:l,k) = sPrev(:,k+1:l)'*yPrev(:,k);
  end
  D = diag(sum(sPrev.*yPrev));
  
  Qty1 = et*(sPrev'*sPrev);
  Qty2 = [et*sPrev, yPrev];
  Qty3 = [Qty1,  L;
          L'  , -D ];
  Hx   = @(x) et*x - Qty2*(Qty3\(Qty2'*x));
  