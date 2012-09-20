function [x, f, output] = pnopt(smoothF, nonsmoothF, x, options)
% pnopt : Proximal Newton-type methods
% 
% [x, f, output] = pnopt(smoothF, nonsmoothF, x) starts at x and seeks a 
%   minimizer of the objective function in composite form. smoothF is a handle
%   to a function that returns the smooth function value and gradient. nonsmoothF
%   is a handle to a function that returns the nonsmooth function value and 
%   proximal mapping. 
% 
% [x, f, output] = pnopt(smoothF, nonsmoothF, x, options) replaces the default
%   optimization parameters with those in options, a structure created using the
%   pnopt_optimset function.
% 
  
% ============ Process options ============
  
  PNoptions = pnopt_optimset(...
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
    options = pnopt_optimset(PNoptions, options);
  else
    options = pnopt_options;
  end
  
  method = options.method;
    
  % ============ Call solver ============
  
  switch method
    case 'Bfgs'
      [x, f, output] = pnopt_PQN(smoothF, nonsmoothF, x, options);
    case 'Lbfgs'
      [x, f, output] = pnopt_PQN(smoothF, nonsmoothF, x, options);
    case 'Newton'
      [x, f, output] = pnopt_PN(smoothF, nonsmoothF, x, options);
    otherwise
      error(sprintf('Unrecognized method ''%s''.', method)) %#ok<SPERR>
  end
  