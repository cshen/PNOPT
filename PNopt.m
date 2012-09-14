function [x, f, output] = ProxQuasiNewton(smoothF, nonsmoothF, x, options)
% ProxQuasiNewton : Proximal quasi-Newton-type methods
% 
% [x, f, output] = ProxQuasiNewton(smoothF, nonsmoothF, x) starts at x and seeks 
%   a minimizer of the objective function in composite form. smoothF is a handle
%   to a function that returns the smooth function value and gradient. nonsmoothF
%   is a handle to a function that returns the nonsmooth function value and 
%   proximal mapping. 
% 
% [x, f, output] = ProxQuasiNewton(smoothF, nonsmoothF, x, options) replaces the   
%   default optimization options with those in options, a structure created 
%   using the PNoptimset function.
% 
  REVISION = '$Revision: 0.1.1$';
  DATE     = '$Date: Sep. 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Process options ============
  
  PNoptions = PNoptimset(...
    'debug'       , 0       ,... % debug mode 
    'descParam'   , 0.0001  ,... % sufficient descent parameter
    'display'     , 10      ,... % display frequency (<= 0 for no display) 
    'LbfgsMem'    , 50      ,... % L-BFGS memory
    'maxfunEvals' , 5000    ,... % max number of function evaluations
    'maxIter'     , 500     ,... % max number of iterations
    'method'      , 'Lbfgs' ,... % method for building Hessian approximation
    'subMethod'   , 'Tfocs' ,... % solver for solving subproblems
    'funTol'      , 1e-9    ,... % stopping tolerance on relative change in the objective function 
    'optTol'      , 1e-6    ,... % stopping tolerance on optimality condition
    'xtol'        , 1e-9     ... % stopping tolerance on solution
    );
  
  if nargin > 3
    options = PNoptimset(PNoptions, options);
  else
    options = PNoptions;
  end
  
  method = options.method;
    
  % ============ Call solver ============
  
  switch method
    case 'Bfgs'
      [x, f, output] = ProxQuasiNewton(smoothF, nonsmoothF, x, options);
    case 'Lbfgs'
      [x, f, output] = ProxQuasiNewton(smoothF, nonsmoothF, x, options);
    case 'Newton'
      [x, f, output] = ProxNewton(smoothF, nonsmoothF, x, options);
    otherwise
      error(sprintf('Unrecognized method ''%s''.', method)) %#ok<SPERR>
  end
  