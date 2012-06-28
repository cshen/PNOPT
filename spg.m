function [x, f, output] = spg(smoothF, nonsmoothF, x, varargin)
% Spg : Spectral proximal gradient methods
% 
% [x, f, output] = Spg(smoothF, nonsmoothF, x) starts at x and seeks a minimizer
%   of the objective function. smoothF is a handle to a function that returns the
%   smooth function value and gradient. nonsmoothF is a handle to a function that
%   returns the nonsmooth function value and prox.
% 
% [x, f, output] = Spg(smoothF, nonsmoothF, x, options) replaces the default 
%   optimization options replaced with values in options, a structure created 
%   using the SetPNoptOptions function.
% 
  REVISION = '$Revision: 0.1.6$';
  DATE     = '$Date: June 24, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ------------ Initialize ------------
  
  n = length(x);

  % Set default options
  defaultOptions = SetPNoptOptions(...
    'checkOpt'        , 1        ,... % Check optimality (requires prox evaluation)
    'display'         , 1        ,... % display level 
    'LSmemory'   , 10       ,... % Number of previous function values to save
    'maxfunEvals'     , 50000    ,... % Max number of function evaluations
    'maxIter'         , 5000     ,... % Max number of iterations
    'printEvery'      , 100      ,... % display output every printEvery iterations
    'TolFun'          , 1e-9     ,... % Stopping tolerance on objective function 
    'TolOpt'          , 1e-6     ,... % Stopping tolerance on optimality
    'TolX'            , 1e-9      ... % Stopping tolerance on solution
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
  if nargin > 3
    options = SetSpgOptions(defaultOptions, varargin{1});
  else
    options = defaultOptions;
  end
  
  checkOpt         = options.checkOpt;
  display          = options.display;
  LSmemory = options.LSmemory;
  maxfunEvals      = options.maxfunEvals;
  maxIter          = options.maxIter;
  printEvery       = options.printEvery;
  TolFun           = options.TolFun;
  TolOpt           = options.TolOpt;
  TolX             = options.TolX;
  
  iter            = 0; 
  Trace.f         = zeros(maxIter+1,1);
  Trace.funEvals  = zeros(maxIter+1,1);
  Trace.proxEvals = zeros(maxIter+1,1);
  if checkOpt
    Trace.optimality = zeros(maxIter+1,1);
  end
  
  if display    
    if checkOpt
      fprintf(' %s\n',repmat('=',1,69));
      fprintf('                     SPG  v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,69));
      fprintf(' %4s  %8s  %10s  %12s  %12s  %12s \n',...
        '','F evals', 'Prox evals', 'Step len.', 'Obj. val.', 'optimality');
      fprintf(' %s\n',repmat('-',1,69));
    else
      fprintf(' %s\n',repmat('=',1,55));
      fprintf('              Spg  v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,55));
      fprintf(' %4s  %8s  %10s  %12s  %12s \n',...
        '','F evals', 'Prox evals', 'Step len.', 'Obj. val.');
      fprintf(' %s\n',repmat('-',1,55));
    end
  end
  
  % ------------ Evaluate objective function at starting x ------------ 
  
  [f, Df] = smoothF(x);
   h      = nonsmoothF(x);
   f      = f + h;
  
  % ------------ Start collecting data for display and output ------------ 
  
  funEvals  = 1;
  proxEvals = 0;
  if checkOpt
    [h1, x1]  = nonsmoothF(x-Df,1); %#ok<ASGLU>
    proxEvals = proxEvals + 1;
    opt       = norm(x1-x,'inf'); 
  end
  
  Trace.f(1)         = f;
  Trace.funEvals(1)  = funEvals;
  Trace.proxEvals(1) = proxEvals;
  if checkOpt
    Trace.optimality(1) = opt; 
  end
  
  if display    
    if checkOpt
      fprintf(' %4d  %8d  %10d  %12s  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, '', f, opt);
    else
      fprintf(' %4d  %8d  %10d  %12s  %12.4e \n',...
        iter, funEvals, proxEvals, '', f);
    end
  end
  
  % ------------ Check if starting x is optimal ------------ 
  
  if checkOpt && opt <= TolOpt
    output = struct(...
      'flag'      , FLAG_OPTIMAL,...
      'funEvals'  , funEvals    ,...
      'Iterations', iter        ,...
      'method'    , method      ,...
      'optimality', opt         ,...
      'options'   , options     ,...
      'proxEvals' , proxEvals   ,...
      'Trace'     , Trace        ...
      );
  
    if display
      fprintf(' %s\n',repmat('-',1,64));
      fprintf(' %s\n',MESSAGE_OPTIMAL)
      fprintf(' %s\n',repmat('-',1,64));
    end
    
    return
  end

% ------------ Main Loop -------------
  
  while 1
    iter = iter+1; 
    
    % ------------ Compute search direction ------------
    
    if iter > 1
      s  = x-xPrev;
      y  = Df-DfPrev;
      BbStepLen  = (y'*s)/(y'*y);
      if BbStepLen <= 1e-9 || 1e9 <= BbStepLen
        BbStepLen = min(1,1/norm(Df,1));
      end
    else
      BbStepLen = min(1,1/norm(Df,1));
    end
    
    % ------------ Conduct line search ------------    
    xPrev   = x;
    if iter+1 > LSmemory
      fPrev = [fPrev(2:end), f];
    else
      fPrev(iter) = f;
    end
    DfPrev  = Df;
    
    % Conduct line search for a step length that safisfies the Armijo condition
    [x, f, Df, step, LSflag ,LSiter] = ...
      CurvySearch(x, -Df, BbStepLen, fPrev, -norm(Df)^2, smoothF, nonsmoothF,...
      TolX, maxfunEvals - funEvals); %#ok<ASGLU>
    
    % ------------ Collect data and display status ------------
    
    funEvals  = funEvals + LSiter;
    proxEvals = proxEvals + LSiter;
    if checkOpt
      [h1, x1] = nonsmoothF(x-Df,1); %#ok<ASGLU>
      proxEvals = proxEvals + 1;
      opt      = norm(x1-x,'inf'); 
    end
    
    Trace.f(iter+1)         = f;
    Trace.funEvals(iter+1)  = funEvals;
    Trace.proxEvals(iter+1) = proxEvals;
    if checkOpt
      Trace.optimality(iter+1) = opt; 
    end
    
    if display && mod(iter,printEvery) == 0
      if checkOpt
        fprintf(' %4d  %8d  %10d  %12.4e  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f, opt);
      else
        fprintf(' %4d  %8d  %10d  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f);
      end
    end
    
    % ------------Check stopping criteria------------
    
    % Check optimality condition
    if checkOpt && opt <= TolOpt
      flag    = FLAG_OPTIMAL;
      Message = MESSAGE_OPTIMAL;
      break
      
    % Check lack of progress
    elseif norm(x-xPrev)/max(1,norm(xPrev)) <= TolX 
      flag    = FLAG_TOLX;
      Message = MESSAGE_TOLX;
      break
    elseif f <= min(fPrev) && abs(min(fPrev)-f)/max(1,abs(fPrev(end))) <= TolFun
      flag    = FLAG_TOLFUN;
      Message = MESSAGE_TOLFUN;
      break
      
    % Check function evaluation/iteration cap
    elseif iter >= maxIter 
      flag    = FLAG_MAXITER;
      Message = MESSAGE_MAXITER;
      break
    elseif funEvals >= maxfunEvals
      flag    = FLAG_MAXFUNEVALS;
      Message = MESSAGE_MAXFUNEVALS;
      break
    end
  end
  
  % -------------------- Cleanup and exit --------------------
  
  Trace.f         = Trace.f(1:iter+1);
  Trace.funEvals  = Trace.funEvals(1:iter+1);
  Trace.proxEvals = Trace.proxEvals(1:iter+1);
  if checkOpt
    Trace.optimality = Trace.optimality(1:iter+1);
  end
  
  if display && mod(iter,printEvery) > 0
    if checkOpt
      fprintf(' %4d  %8d  %10d  %12.4e  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f, opt);
    else
      fprintf(' %4d  %8d  %10d  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f);
    end
  end
  
  output = struct(...
    'flag'      , flag     ,...
    'funEvals'  , funEvals ,...
    'Iterations', iter     ,...
    'options'   , options  ,...
    'proxEvals' , proxEvals,...
    'Trace'     , Trace     ...
    );
  if checkOpt
     output.Optimality = opt;
  end
  
  if display
    if checkOpt
      fprintf(' %s\n',repmat('-',1,57));
      fprintf(' %s\n',Message)
      fprintf(' %s\n',repmat('-',1,57));
    else
      fprintf(' %s\n',repmat('-',1,43));
      fprintf(' %s\n',Message)
      fprintf(' %s\n',repmat('-',1,43));
    end
  end
