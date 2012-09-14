function [f, Df] = smooth_quad(Hf, Df, f, x)
% smooth_quad : Quadratic function
%
%   $Revision: 0.1.2 $  $Date: 2012/06/24 $
% 
  global eff quadf quadDf
  
  if isa(Hf,'function_handle')
    Hfx = Hf(x);
  elseif isnumeric(Hf)
    Hfx = Hf*x;
  end
  
  f = 0.5*x'*Hfx + Df'*x + f;
  if nargout > 1
    quadDf = Hfx + Df;
    Df = quadDf;
  end