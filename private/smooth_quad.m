function [ f_y, Df_y ] = smooth_quad( P, q, r, x )
% smooth_quad : Quadratic function
%
%   $Revision: 0.1.2 $  $Date: 2012/09/15 $
% 
  global subProb_Dg_y
  
  if isa( P, 'function_handle' )
    H_x = P( x );
  else
    H_x = P * x;
  end
  
  f_y = 0.5 * x' * H_x + q' * x + r;
  if nargout > 1
    subProb_Dg_y = H_x + q;
    Df_y = subProb_Dg_y;
  end