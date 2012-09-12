function [f, Df] = pnopt_subDual(d, B, x, u, nonsmoothF, conjF)
% pnopt_subDual : Dual of the proximal Newton subproblem
% 
%   $Revision: 0.1.2 $  $Date: 2012/07/15 $
% 
% ============ unscaled dual ============
% 
%   [~, v] = nonsmoothF(x-u/d,1/d);
%   if isa(B, 'function_handle')
%     Bu = B(u);
%     Df = -v + Bu + x;
%   elseif isnumeric(B)
%     Bu = B\u;
%     Df = -v + Bu + x;
%   end
%   
%   [fv, v] = conjF(d*x-u,d);
%        f  = fv + 0.5*u'*Bu + x'*u + 0.5/d*norm(d*x-u-v)^2;
%   
% ============ scaled dual ============
  
  d2 = sqrt(d);
  
  conjF_scale      = prox_scale( conjF, d2 );
  nonsmoothF_scale = prox_scale( nonsmoothF, 1/d2 );
  
  [~, v] = nonsmoothF_scale(d2*x-u,1);
  if isa(B, 'function_handle')
    Bu = B(u);
    Df = -v + d*Bu + d2*x;
  elseif isnumeric(B)
    Bu = B\u;
    Df = -v + d*Bu + d2*x;
  end
  
  [fv, v] = conjF_scale(d2*x-u,1);
       f  = fv + 0.5*norm(d2*x-u-v)^2 + 0.5*d*u'*Bu + d2*x'*u;
  