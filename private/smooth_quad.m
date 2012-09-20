function [ f_x, Df_x ] = smooth_quad( Hf_x, Df_x, f_x, x )
% smooth_quad : Quadratic function
%
%   $Revision: 0.1.2 $  $Date: 2012/06/24 $
% 
  global quad_Df_x
  
  if isa( Hf_x, 'function_handle' )
    Hfx = Hf_x( x );
  elseif isnumeric( Hf_x )
    Hfx = Hf_x * x;
  end
  
  f_x = 0.5 * x' * Hfx + Df_x' * x + f_x;
  if nargout > 1
    quad_Df_x = Hfx + Df_x;
    Df_x = quad_Df_x;
  end