function [M, Dh] = smoothDual(d, Q, x, nonsmoothF, v)
% smoothDual : Smoothed dual
% 
%   $Revision: 0.1.0 $  $Date: 2012/06/24 $
% 
  [~, w] = nonsmoothF(x-v/d,1/d);
  if isa(Q, 'function_handle')
    Dh   = -w + Q(v) + x;
  elseif nnz(tril(Q,-1)) == 0
    Dh   = -w + Q\(Q'\v) + x;
  else
    error('smoothDual:BadQuad', 'Second argument must be a function handle or Cholesky factor')
  end
  M      = 0.5*norm(Dh)^2;
  