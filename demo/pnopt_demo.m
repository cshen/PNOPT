% pnopt_demo : PNOPT demo using the Graphical Lasso on flow cytometry data
% 
%   $Revision: 0.0.1 $  $Date: 2012/11/10 $
% 
%% Load flow cytometry data
  
  load flow.mat
  
  C = 
  
  log_likelihood = @(T) smooth_logdet( 0.5, inv(X) );
  
  lambda    = 15*ones(n+1,1);
  lambda(1) = 0;      % Do not penalize bias term
  L1pen     = prox_l1(lambda);
  
  w0 = zeros(p,1);
  
%% Solve the Graphical Lasso using PNOPT
  
  options = PNoptimset('debug', 1);
  [w, f, output] = pnopt(log_likelihood, L1pen, w0, options);
  