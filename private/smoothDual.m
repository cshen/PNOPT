function [M, Dh] = smoothedF(d, Be, x, nonsmoothF, v)
% smoothedF : Smoothed (Moreau-Yosida) 
%
%   $Revision: 0.1.0 $  $Date: 2012/06/24 $
% 
  [~, w] = nonsmoothF(x-v/d,1/d);
  if isa(Be, 'function_handle')
    Dh   = -w + Be(v) + x;
  elseif nnz(tril(Be,-1)) == 0
    Dh   = -w + Be\(Be'\v) + x;
  else
    error('PNdual:BadBe', 'Second argument must be a function handle or Cholesky factor')
  end
  M      = 0.5*norm(Dh)^2;