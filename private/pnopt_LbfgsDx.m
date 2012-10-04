function Dx = pnopt_LbfgsDx( s_old, y_old, de, Df_x ) %#codegen
% pnopt_LbfgsDx : Product with L-BFGS Hessian approximation inverse
% 
%   $Revision: 0.1.1 $  $Date: 2012/06/30 $
% 
  qty1 = 1 ./ sum( y_old.*s_old );
  l = length( qty1 );
  
  % two loop recursion
  qty2 = zeros( l, 1 );
  Dx   = - Df_x;
  for k = l:-1:1
    qty2(k) = qty1(k) * s_old(:,k)' * Dx;
    Dx      = Dx - qty2(k) * y_old(:,k);
  end
  qty3 = zeros(l,1);
  Dx   = Dx / de;
  for k = 1:l
    qty3(k) = qty1(k) * y_old(:,k)' * Dx;
    Dx      = Dx + s_old(:,k) * ( qty2(k) - qty3(k) );
  end
  