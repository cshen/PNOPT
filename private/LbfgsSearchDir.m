function Dx = LbfgsSearchDir(sPrev, yPrev, et, Df) %#codegen
% LbfgsSearchDir : L-BFGS search direction
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $
% 
  qty1 = 1./sum(yPrev.*sPrev);
  l = length(qty1);
  
  % two loop recursion
  qty2 = zeros(l,1);
  Dx = -Df;
  for k = l:-1:1
    qty2(k) = qty1(k)*sPrev(:,k)'*Dx;
    Dx = Dx - qty2(k)*yPrev(:,k);
  end
  qty3 = zeros(l,1);
  Dx = Dx/et;
  for k = 1:l
    qty3(k) = qty1(k)*yPrev(:,k)'*Dx;
    Dx = Dx + sPrev(:,k)*(qty2(k) - qty3(k));
  end
  