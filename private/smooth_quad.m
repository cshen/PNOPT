function [f, Df] = smooth_quad(P, q, x)
% smooth_quad : Quadratic function
%
%   $Revision: 0.1.2 $  $Date: 2012/06/24 $
% 
  global quadDf
  
  if isa(P,'function_handle')
    Px = P(x);
  elseif isnumeric(P)
    Px = P*x;
  end
  
  f = 0.5*x'*Px + q'*x;
  if nargout > 1
    quadDf = Px + q;
    Df = quadDf;
  end
  