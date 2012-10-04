function [ f_y, Df_y ] = pnopt_quad( P, q, r, x )
% pnopt_quad : Generalized proximal term
%
%   $Revision: 0.1.2 $  $Date: 2012/06/24 $
% 
  global subprob_Df_y
  
  if isa( P, 'function_handle' )
    H_x = P( x );
  elseif isnumeric( P )
    H_x = P * x;
  end
  
  f_y = 0.5 * x' * H_x + q' * x + r;
  if nargout > 1
    subprob_Df_y = H_x + q;
    Df_y = subprob_Df_y;
  end