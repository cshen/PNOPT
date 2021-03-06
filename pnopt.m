function [ x, f, output ] = pnopt( smoothF, nonsmoothF, x, options )
% pnopt : Proximal Newton-type methods
% 
% [ x, f, output ] = pnopt( smoothF, nonsmoothF, x ) starts at x and seeks a 
%   minimizer of the objective function in composite form. smoothF is a handle
%   to a function that returns the smooth function value and gradient. nonsmoothF
%   is a handle to a function that returns the nonsmooth function value and 
%   proximal mapping. 
%  
% [ x, f, output ] = pnopt( smoothF, nonsmoothF, x, options ) replaces the default
%   optimization parameters with those in options, a structure created using the
%   pnopt_optimset function.
% 
%   $Revision: 0.8.0 $  $Date: 2012/12/01 $
  
% ============ Process options ============
  
  if exist( 'tfocs', 'file')
    default_options = pnopt_optimset( ...
      'debug'          , 0       ,... % debug mode 
      'desc_param'     , 0.0001  ,... % sufficient descent parameter
      'display'        , 10      ,... % display frequency (<= 0 for no display) 
      'Lbfgs_mem'      , 50      ,... % L-BFGS memory
      'maxfunEv'       , 5000    ,... % max number of function evaluations
      'maxIter'        , 500     ,... % max number of iterations
      'method'         , 'Lbfgs' ,... % method for building Hessian approximation
      'subProb_solver' , 'tfocs' ,... % solver for solving subproblems
      'ftol'           , 1e-9    ,... % stopping tolerance on relative change in the objective function 
      'optim_tol'      , 1e-6    ,... % stopping tolerance on optimality condition
      'xtol'           , 1e-9     ... % stopping tolerance on solution
      );
  else
    default_options = pnopt_optimset( ...
      'debug'          , 0        ,... 
      'desc_param'     , 0.0001   ,... 
      'display'        , 10       ,... 
      'Lbfgs_mem'      , 50       ,... 
      'maxfunEv'       , 5000     ,... 
      'maxIter'        , 500      ,... 
      'method'         , 'Lbfgs'  ,... 
      'subProb_solver' , 'sparsa' ,... 
      'ftol'           , 1e-9     ,... 
      'optim_tol'      , 1e-6     ,...
      'xtol'           , 1e-9      ... 
      );
  end
  
  if nargin > 3
    options = pnopt_optimset( default_options, options );
  else
    options = pnopt_options;
  end
  
  method = options.method;
    
  % ============ Call solver ============
  
  switch method
    case {'bfgs', 'Lbfgs'}
      [ x, f, output ] = pnopt_PQN( smoothF, nonsmoothF, x, options );
    case 'newton'
      [ x, f, output ] = pnopt_PN( smoothF, nonsmoothF, x, options );
    otherwise
      error( sprintf( 'Unrecognized method ''%s''.', method ) ) %#ok<SPERR>
  end
  