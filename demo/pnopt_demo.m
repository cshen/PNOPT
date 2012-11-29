% pnopt_demo : PNOPT demo using the Graphical Lasso on flow cytometry data
% 
%   $Revision: 0.0.1 $  $Date: 2012/11/10 $
% 
%% Load flow cytometry data
  
  load flow.mat
  
  loglikelihood = @(T) smooth_logdet();
  
%% Solve the Graphical Lasso using PNOPT
  
  lambda    = 15*ones(n+1,1);
  lambda(1) = 0;      % Do not penalize bias term
  L1pen     = prox_l1(lambda);
  
  w0 = zeros(p,1);
  
%% 
  
  PNoptions = PNoptimset('debug', 1);
  [wPN, fPN, PNoutput] = ProxQuasiNewton(loglikelihood, L1pen, w0, PNoptions);
  