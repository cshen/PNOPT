function [f, Df] = smooth_quad(Hx, g, x, y)
% smooth_quad : Quadratic function
%
%   $Revision: 0.1.2 $  $Date: 2012/06/24 $
% 
  global quadDf
  
  dx  = y-x;
  if isa(Hx,'function_handle')
    Hdx = Hx(dx);
  elseif isnumeric(Hx)
    Hdx = Hx*dx;
  end
  
  f = 0.5*dx'*Hdx + g'*dx;
  if nargout > 1
    quadDf = Hdx + g;
    Df = quadDf;
  end
  