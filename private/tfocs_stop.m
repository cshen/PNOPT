function stop = tfocs_stop( x, nonsmoothF, optTol ) 
% tfocs_stop : TFOCS stopping condition
%
%   $Revision: 0.1.2 $  $Date: 2012/09/15 $
% 
  global subProb_Dg_y subProb_optim
  
  [ ~, x_prox ]   = nonsmoothF( x - subProb_Dg_y ,1);
    subProb_optim = norm( x_prox - x ,'inf');
    stop          = subProb_optim <= optTol;
  
    