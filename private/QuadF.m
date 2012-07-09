function [f, Df] = QuadF(Hx, g, x, y)
% QuadF : Quadratic function
%
%   $Revision: 0.1.2 $  $Date: 2012/06/24 $

  dx  = y-x;
  if isa(Hx,'function_handle')
    Hdx = Hx(dx);
  elseif isnumeric(Hx)
    Hdx = Hx*dx;
  else
    error('QuadF:badHess', 'Hessian must be a function handle or numeric array.')
  end
  f   = 0.5*dx'*Hdx + g'*dx;
  Df  = Hdx + g;
  
  