function stop = tfocs_stop(x, nonsmoothF, tolopt) 
% tfocs_stop : TFOCS stopping condition
%
%   $Revision: 0.1.0 $  $Date: 2012/07/10 $
% 
  global quadDf quadopt
  
  [~, xProx ] = nonsmoothF( x - quadDf ,1);
     quadopt  = norm( xProx - x ,'inf');
     stop     = quadopt <= tolopt;
  