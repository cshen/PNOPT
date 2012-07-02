function [x, f, output] = spg(smoothF, nonsmoothF, x, varargin)
% Spg : Spectral proximal gradient methods
% 
% [x, f, output] = Spg(smoothF, nonsmoothF, x) starts at x and seeks a minimizer
%   of the objective function in composite form. smoothF is a handle to a function
%   that returns the smooth function value and gradient. nonsmoothF is a handle 
%   to a function that returns the nonsmooth function value and prox.
% 
% [x, f, output] = Spg(smoothF, nonsmoothF, x, options) replaces the default 
%   optimization options with those in options, a structure created using the 
%   SetPNoptOptions function.
% 
  REVISION = '$Revision: 0.1.6$';
  DATE     = '$Date: June 24, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Initialize ============
  
  % Set default options
  defaultOptions = SetPNoptOptions(...
    'checkOpt'    , 1     ,... % Check optimality (requires prox evaluation)
    'display'     , 1     ,... % display frequency (<= 0 for no display) 
    'LSmemory'    , 10    ,... % Number of previous function values to save
    'maxfunEvals' , 50000 ,... % Max number of function evaluations
    'maxIter'     , 5000  ,... % Max number of iterations
    'TolFun'      , 1e-9  ,... % Stopping tolerance on objective function 
    'TolOpt'      , 1e-6  ,... % Stopping tolerance on optimality
    'TolX'        , 1e-9   ... % Stopping tolerance on solution
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
  
  checkOpt    = options.checkOpt;
  display     = options.display;
  LSmemory    = options.LSmemory;
  maxfunEvals = options.maxfunEvals;
  maxIter     = options.maxIter;
  TolFun      = options.TolFun;
  TolOpt      = options.TolOpt;
  TolX        = options.TolX;
  
  iter            = 0; 
  loop            = 1; 
  Trace.f         = zeros(maxIter+1,1);
  Trace.funEvals  = zeros(maxIter+1,1);
  Trace.proxEvals = zeros(maxIter+1,1);
  if checkOpt
    Trace.optimality = zeros(maxIter+1,1);
  end
  
  if display > 0    
    if checkOpt
      fprintf(' %s\n',repmat('=',1,64));
      fprintf('            SPG  v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,64));
      fprintf(' %4s   %6s  %6s  %12s  %12s  %12s \n',...
        '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'optimality');
      fprintf(' %s\n',repmat('-',1,64));
    else
      fprintf(' %s\n',repmat('=',1,50));
      fprintf('     SPG  v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,50));
      fprintf(' %4s   %6s  %6s  %12s  %12s \n',...
        '','Fun.', 'Prox', 'Step len.', 'Obj. val.');
      fprintf(' %s\n',repmat('-',1,50));
    end
  end
  
  % ============ Evaluate objective function at starting x ============ 
  
  [f, Df] = smoothF(x);
   h      = nonsmoothF(x);
   f      = f + h;
  
  % ============ Start collecting data for display and output ============ 
  
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
  
  if display > 0    
    if checkOpt
      fprintf(' %4d | %6d  %6d  %12s  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, '', f, opt);
    else
      fprintf(' %4d | %6d  %6d  %12s  %12.4e \n',...
        iter, funEvals, proxEvals, '', f);
    end
  end
  
  % ============ Check if starting x is optimal ============ 
  
  if checkOpt && opt <= TolOpt
    flag    = FLAG_OPTIMAL;
    message = MESSAGE_OPTIMAL;
    loop    = 0;
  end

% ============ Main Loop ============
  
  while loop
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
    
    if display > 0 && mod(iter,display > 0) == 0
      if checkOpt
        fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f, opt);
      else
        fprintf(' %4d | %6d  %6d  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f);
      end
    end
    
    % ------------ Check stopping criteria ------------
    
    % Check optimality condition
    if checkOpt && opt <= TolOpt
      flag    = FLAG_OPTIMAL;
      message = MESSAGE_OPTIMAL;
      loop    = 0;
      
    % Check lack of progress
    elseif norm(x-xPrev,'inf')/max(1,norm(xPrev,'inf')) <= TolX 
      flag    = FLAG_TOLX;
      message = MESSAGE_TOLX;
      loop    = 0;
    elseif f <= min(fPrev) && abs(min(fPrev)-f)/max(1,abs(fPrev(end))) <= TolFun
      flag    = FLAG_TOLFUN;
      message = MESSAGE_TOLFUN;
      loop    = 0;
      
    % Check function evaluation/iteration cap
    elseif iter >= maxIter 
      flag    = FLAG_MAXITER;
      message = MESSAGE_MAXITER;
      loop    = 0;
    elseif funEvals >= maxfunEvals
      flag    = FLAG_MAXFUNEVALS;
      message = MESSAGE_MAXFUNEVALS;
      loop    = 0;
    end
  end
  
  % ============ Cleanup and exit ============
  
  Trace.f         = Trace.f(1:iter+1);
  Trace.funEvals  = Trace.funEvals(1:iter+1);
  Trace.proxEvals = Trace.proxEvals(1:iter+1);
  if checkOpt
    Trace.optimality = Trace.optimality(1:iter+1);
  end
  
  if display > 0 && mod(iter,display > 0) > 0
    if checkOpt
      fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f, opt);
    else
      fprintf(' %4d | %6d  %6d  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f);
    end
  end
  
  output = struct(...
    'flag'      , flag      ,...
    'funEvals'  , funEvals  ,...
    'iterations', iter      ,...
    'options'   , options   ,...
    'proxEvals' , proxEvals ,...
    'Trace'     , Trace      ...
    );
  if checkOpt
     output.Optimality = opt;
  end
  
  if display > 0
    if checkOpt
      fprintf(' %s\n',repmat('-',1,64));
      fprintf(' %s\n',message)
      fprintf(' %s\n',repmat('-',1,64));
    else
      fprintf(' %s\n',repmat('-',1,50));
      fprintf(' %s\n',message)
      fprintf(' %s\n',repmat('-',1,50));
    end
  end
