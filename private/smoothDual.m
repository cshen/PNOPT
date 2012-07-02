function [M, Dh] = smoothDual(d, B, x, nonsmoothF, v)
% smoothDual : Smoothed dual
% 
%   $Revision: 0.1.0 $  $Date: 2012/06/24 $
% 
  [~, w] = nonsmoothF(x-v/d,1/d);
  if isa(B, 'function_handle')
    Dh   = -w + B(v) + x;
  elseif nnz(tril(Be,-1)) == 0
    Dh   = -w + B\(B'\v) + x;
  else
    error('smoothDual:BadB', 'Second argument must be a function handle or Cholesky factor')
  end
  M      = 0.5*norm(Dh)^2;
  