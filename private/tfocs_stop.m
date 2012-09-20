function stop = tfocs_stop( x, nonsmoothF, optTol ) 
% tfocs_stop : TFOCS stopping condition
%
%   $Revision: 0.1.0 $  $Date: 2012/07/10 $
% 
  global quad_Df_x quad_opt
  
  [ ~, x_prox ] = nonsmoothF( x - quad_Df_x ,1);
    quad_opt    = norm( x_prox - x ,'inf');
    stop        = quad_opt <= optTol;
  