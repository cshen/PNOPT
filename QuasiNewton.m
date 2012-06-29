function [x, f, output] = QuasiNewton(Fun, x, varargin)
% QN : Quasi-Newton Methods
% 
% [x, f, output] = QuasiNewton(fun, x) starts at x and seeks a minimizer of the
%   objective function. fun is a handle to a function that returns the objective
%   function value and gradient.
% 
% [x, f, output] = QuasiNewton(fun, x, options) replaces the default 
%   optimization options replaced with values in options, a struct created using
%   the SetNewtonOptions function.
% 
  REVISION = '$Revision: 0.2.0$';
  DATE     = '$Date: June 24, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ------------ Initialize ------------
  
  % Set default options
  defaultOptions = SetPNoptOptions(...
    'display'         , 1      ,... % display level
    'LbfgsCorrections', 50     ,... % Number of L-BFGS corrections
    'maxfunEvals'     , 50000  ,... % Max number of function evaluations
    'maxIter'         , 500    ,... % Max number of iterations
    'method'          , 'Lbfgs',... % method for choosing search directions
    'printEvery'      , 10     ,... % display output every printEvery iterations
    'TolFun'          , 1e-9   ,... % Stopping tolerance on objective function 
    'TolOpt'          , 1e-6   ,... % Stopping tolerance on optimality
    'TolX'            , 1e-9    ... % Stopping tolerance on solution
    );
  
  % Set stopping flags and messages
  FLAG_OPTIMAL     = 1;
  FLAG_TOLX        = 2;
  FLAG_TOLFUN      = 3;
  FLAG_MAXITER     = 4;
  FLAG_MAXFUNEVALS = 5;
  
  MESSAGE_OPTIMAL     = 'Optimality below TolOpt.';
  MESSAGE_TOLX        = 'Relative change in x below TolX.';
  MESSAGE_TOLFUN      = 'Relative change in function value below TolFun.';
  MESSAGE_MAXITER     = 'Max number of iterations reached.';
  MESSAGE_MAXFUNEVALS = 'Max number of function evaluations reached.';
  
  % Replace default option values with values in user-supplied options struct
  if nargin > 2
    options = SetPNoptOptions(defaultOptions, varargin{1});
  else
    options = defaultOptions;
  end
  
  LbfgsCorrections = options.LbfgsCorrections;
  display          = options.display;
  maxfunEvals      = options.maxfunEvals;
  maxIter          = options.maxIter;
  method           = options.method;
  printEvery       = options.printEvery;
  TolFun           = options.TolFun;
  TolOpt           = options.TolOpt;
  TolX             = options.TolX;
    
  if strcmp(method,'Newton')
    Hess = options.Hess;
  end
  
  iter             = 0; 
  Trace.f          = zeros(maxIter+1,1);
  Trace.funEvals   = zeros(maxIter+1,1);
  Trace.optimality = zeros(maxIter+1,1);
  
  if display
    fprintf(' %s\n',repmat('=',1,56));
    fprintf('          QuasiNewton  v.%s (%s)\n', REVISION, DATE);
    fprintf(' %s\n',repmat('=',1,56));
    fprintf(' %4s   %6s  %12s  %12s  %12s \n',...
      '','Fun.', 'Step len.', 'Obj. val.', 'optimality');
    fprintf(' %s\n',repmat('-',1,56));
  end
  
  % ------------ Evaluate objective function at starting x ------------ 
  
  [f, Df] = Fun(x);
  if strcmp(method,'Newton')  % Evaluate Hessian if necessary.
    Hf    = Hess(x);
    if ~isa(Hf,'numeric') && ~isa(Hf,'function_handle')
      error('QuasiNewton:BadHess', 'options.Hess must be a function handle or numeric array')
    end
  end
  
  % ------------ Start collecting data for display and output ------------ 
  
  funEvals = 1;
  opt      = norm(Df,'inf');
  
  Trace.f(1)          = f;
  Trace.funEvals(1)   = funEvals;
  Trace.optimality(1) = opt;
  
  if display
    fprintf(' %4d | %6d  %12s  %12.4e  %12.4e\n',...
      iter, funEvals, '', f, opt);
  end
  
  % ------------ Check if starting x is optimal ------------ 
  
  if opt <= TolOpt      
    output = struct(...
    'flag'      , FLAG_OPTIMAL,...
    'funEvals'  , funEvals    ,...
    'iterations', iter        ,...
    'method'    , method      ,...
    'optimality', opt         ,...
    'options'   , options     ,...
    'Trace'     , Trace        ...
      );
  
    if display
      fprintf(' %s\n',repmat('-',1,56));
      fprintf(' %s\n',MESSAGE_OPTIMAL);
      fprintf(' %s\n',repmat('-',1,56));
    end
  
    return
  end
  
% ------------ Main Loop -------------
  
  while 1
    iter = iter+1; 
    
    % ------------ Compute search direction ------------
    
    switch method
      % BFGS method
      case 'bfgs'
        if iter > 1
          s = x-xPrev;
          y = Df-DfPrev;
          qty1 = R'*(R*s);
          if s'*y > 1e-9
            R = cholupdate(cholupdate(R, y/sqrt(y'*s)), qty1/sqrt(s'*qty1),'-');
          end
          Dx = -R\(R'\Df);
        else
          R  = eye(length(x));
          Dx = -Df;
        end
        
      % Limited-memory BFGS method
      case 'Lbfgs'
        if iter > 1
          s = x-xPrev;
          y = Df-DfPrev;
          if y'*s > 1e-9
            if size(sPrev,2) > LbfgsCorrections
              sPrev = [sPrev(:,2:LbfgsCorrections), s];
              yPrev = [yPrev(:,2:LbfgsCorrections), y];
              et    = (y'*y)/(y'*s);
            else
              sPrev = [sPrev, s]; %#ok<AGROW>
              yPrev = [yPrev, y]; %#ok<AGROW>
              et    = (y'*y)/(y'*s);
            end
          end
          
          Dx = LbfgsSearchDir(sPrev, yPrev, et, Df);
        else
          sPrev = zeros(length(x), 0);
          yPrev = zeros(length(x), 0);
          et    = 1;
          Dx    = -Df;
        end
        
      % Newton's method
      case 'Newton'
        % If Hf is a function handle, use pcg to solve Newton system inexactly.
        if isa(Hf,'function_handle')
          Dx = pcg(Hf, -Df, min(0.5,sqrt(opt))*opt);
        elseif isa(Hf,'numeric')
          Dx = Hf\(-Df);
        end
    end
    
    % ------------ Conduct line search ------------
    
    xPrev  = x;
    fPrev  = f;
    DfPrev = Df;
    
    % Conduct line search for a step length that safisfies the Wolfe conditions
    if iter > 1 
      [x, f, Df, step, ~, LSiter] = ...
        cvsrch(Fun, x, f, Df, Dx, 1, max(TolX,1e-9), maxfunEvals - funEvals);
    else
      [x, f, Df, step, ~, LSiter] = ...
        cvsrch(Fun, x, f, Df, Dx, min(1,1/norm(Df,1)), max(TolX,1e-9), maxfunEvals - funEvals);
    end
    
    if strcmp(method,'Newton')  % Evaluate Hessian if necessary.
      Hf = Hess(x);
      if ~isa(Hf,'numeric') && ~isa(Hf,'function_handle') 
        error('QuasiNewton:BadHess', 'options.Hess must be a function handle or numeric array')
      end
    end
    
    % ------------ Collect data for display and output ------------
    
    funEvals   = funEvals + LSiter;   
    opt = norm(Df,'inf');
    
    Trace.f(iter+1)          = f;
    Trace.funEvals(iter+1)   = funEvals;
    Trace.optimality(iter+1) = opt;
    
    if display && mod(iter,printEvery) == 0
      fprintf(' %4d | %6d  %12.4e  %12.4e  %12.4e\n',...
        iter, funEvals, step, f, opt);
    end
    
    % ------------ Check stopping criteria ------------
    
    % Check optimality condition
    if opt <= TolOpt
      flag    = FLAG_OPTIMAL;
      message = MESSAGE_OPTIMAL;
      break
      
    % Check lack of progress
    elseif norm(x-xPrev)/max(1,norm(xPrev)) <= TolX 
      flag    = FLAG_TOLX;
      message = MESSAGE_TOLX;
      break
    elseif f <= fPrev && (fPrev-f)/max(1,abs(fPrev)) <= TolFun
      flag    = FLAG_TOLFUN;
      message = MESSAGE_TOLFUN;
      break
      
    % Check function evaluation/iteration cap
    elseif iter >= maxIter 
      flag    = FLAG_MAXITER;
      message = MESSAGE_MAXITER;
      break
    elseif funEvals >= maxfunEvals
      flag    = FLAG_MAXFUNEVALS;
      message = MESSAGE_MAXFUNEVALS;
      break
    end
  end
  
  % -------------------- Cleanup and exit --------------------
  
  Trace.f          = Trace.f(1:iter+1);
  Trace.funEvals   = Trace.funEvals(1:iter+1);
  Trace.optimality = Trace.optimality(1:iter+1);
  
  if display && mod(iter,printEvery) > 0
    fprintf(' %4d | %6d  %12.4e  %12.4e  %12.4e\n',...
      iter, funEvals, step, f, opt);
  end
      
  output = struct(...
    'flag'      , flag    ,...
    'funEvals'  , funEvals,...
    'iterations', iter    ,...
    'method'    , method  ,...
    'optimality', opt     ,...
    'options'   , options ,...
    'Trace'     , Trace    ...
    );
  
  if display
    fprintf(' %s\n',repmat('-',1,56));
    fprintf(' %s\n',message)
    fprintf(' %s\n',repmat('-',1,56));
  end
