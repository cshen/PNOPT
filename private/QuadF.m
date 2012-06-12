function [f, Df] = QuadF(Hx, g, x, y)
% QuadF : Quadratic function
%
%   $Revision: 0.1.0 $  $Date: 2012/05/01 $

  dx  = y-x;
  Hdx = Hx(dx);
  f   = 0.5*dx'*Hdx + g'*dx;
  Df  = Hdx + g;
  