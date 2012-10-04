function stop = tfocs_stop( x, nonsmoothF, optTol ) 
% tfocs_stop : TFOCS stopping condition
%
%   $Revision: 0.1.1 $  $Date: 2012/07/15 $
% 
  global subprob_Df_y subprob_optim
  
  [ ~, x_prox ]   = nonsmoothF( x - subprob_Df_y ,1);
    subprob_optim = norm( x_prox - x ,'inf');
    stop          = subprob_optim <= optTol;
  