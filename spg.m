function [x, f, output] = spg(smoothF, nonsmoothF, x, varargin)
% spg : Spectral proximal gradient methods
% 
% [x, f, output] = spg(smoothF, nonsmoothF, x) starts at x and seeks a minimizer
%   of the objective function. smoothF is a handle to a function that returns the
%   smooth function value and gradient. nonsmoothF is a handle to a function that
%   returns the nonsmooth function value and prox.
% 
% [x, f, output] = spg(smoothF, nonsmoothF, x, options) replaces the default 
%   optimization options replaced with values in options, a structure created 
%   using the SetPNoptOptions function.
% 
  REVISION = '$Revision: 0.1.2$';
  DATE     = '$Date: June 12, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% -------------------- Initialize --------------------
  
  n = length(x);

  % Set default options
  DefaultOptions = SetPNoptOptions(...
    'CheckOpt'        , 1        ,... % Check optimality (requires prox evaluation)
    'Display'         , 1        ,... % Display level 
    'LineSearchMemory', 10       ,... % Number of previous function values to save
    'MaxFunEvals'     , 50000    ,... % Max number of function evaluations
    'MaxIter'         , 5000     ,... % Max number of iterations
    'PrintEvery'      , 100      ,... % Display output every PrintEvery iterations
    'SmallFirstStep'  , 1        ,... % Choose small step length during first iteration
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
    options = SetSpgOptions(DefaultOptions, varargin{1});
  else
    options = DefaultOptions;
  end
  
  CheckOpt         = options.CheckOpt;
  Display          = options.Display;
  LineSearchMemory = options.LineSearchMemory;
  MaxFunEvals      = options.MaxFunEvals;
  MaxIter          = options.MaxIter;
  PrintEvery       = options.PrintEvery;
  SmallFirstStep   = options.SmallFirstStep;
  TolFun           = options.TolFun;
  TolOpt           = options.TolOpt;
  TolX             = options.TolX;
  
  iter = 0; 
  
  % Evaluate objective function at starting x
  [f, Df] = smoothF(x);
   h      = nonsmoothF(x);
   f      = f + h;
  
  % Start collecting data for display and output
  FunEvals = 1;
  if CheckOpt
    [h1, x1] = nonsmoothF(x-Df,1); %#ok<ASGLU>
    opt      = norm(x1-x)/sqrt(n); 
  end
  
  Trace.f        = zeros(MaxIter+1,1);
  Trace.FunEvals = zeros(MaxIter+1,1);
  if CheckOpt
    Trace.Optimality = zeros(MaxIter+1,1);
  end
  
  Trace.f(1)        = f;
  Trace.FunEvals(1) = FunEvals;
  if CheckOpt
    Trace.Optimality(1) = opt; 
  end
  
  if Display    
    if CheckOpt
      fprintf(' %s\n',repmat('=',1,57));
      fprintf('             SPG  v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,57));
      fprintf(' %4s  %8s  %12s  %12s  %12s \n',...
        '','F evals', 'Step len.', 'Obj. val.', 'Optimality');
      fprintf(' %s\n',repmat('-',1,57));
      fprintf(' %4d  %8d  %12s  %12.4e  %12.4e\n',...
        iter, FunEvals, '', f, opt);
    else
      fprintf(' %s\n',repmat('=',1,43));
      fprintf('       SPG  v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,43));
      fprintf(' %4s  %8s  %12s  %12s  %12s \n',...
        '','F evals', 'Step len.', 'Obj. val.');
      fprintf(' %s\n',repmat('-',1,43));
      fprintf(' %4d  %8d  %12s  %12.4e \n',...
        iter, FunEvals, '', f);
    end
  end
  
  % Check if starting x is optimal
  if CheckOpt && opt <= TolOpt
    output = struct(...
      'Flag'            , FLAG_OPTIMAL    ,...
      'FunEvals'        , FunEvals        ,...
      'Iterations'      , iter            ,...
      'Optimality'      , opt             ,...
      'Trace'           , Trace            ...
      );
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
      if SmallFirstStep
        BbStepLen = min(1,1/norm(Df,1));
      else
        BbStepLen = 1;
      end
    end
    
    % ------------ Conduct line search ------------    
    xPrev   = x;
    if iter+1 > LineSearchMemory
      fPrev = [fPrev(2:end), f];
    else
      fPrev(iter) = f;
    end
    DfPrev  = Df;
    
    % Conduct line search for a step length that safisfies the Armijo condition
    [x, f, Df, step, LineSearchFlag ,LineSearchFunEvals] = ...
      CurvySearch(x, -Df, BbStepLen, fPrev, -norm(Df)^2, smoothF, nonsmoothF,...
      TolX/max(1,norm(Df)), MaxFunEvals-FunEvals); %#ok<ASGLU>
    
    % ------------ Collect data and display status ------------
    
    FunEvals = FunEvals + LineSearchFunEvals;
    if CheckOpt
      [h1, x1] = nonsmoothF(x-Df,1); %#ok<ASGLU>
      opt      = norm(x1-x)/sqrt(n); 
    end
    
    Trace.f(iter+1)        = f;
    Trace.FunEvals(iter+1) = FunEvals;
    if CheckOpt
      Trace.Optimality(iter+1) = opt; 
    end
    
    if Display && mod(iter,PrintEvery) == 0
      if CheckOpt
        fprintf(' %4d  %8d  %12.4e  %12.4e  %12.4e\n',...
          iter, FunEvals, step, f, opt);
      else
        fprintf(' %4d  %8d  %12.4e  %12.4e\n',...
          iter, FunEvals, step, f);
      end
    end
    
    % ------------Check stopping criteria------------
    
    % Check optimality condition
    if CheckOpt && opt <= TolOpt
      Flag    = FLAG_OPTIMAL;
      Message = MESSAGE_OPTIMAL;
      break
      
    % Check lack of progress
    elseif norm(x-xPrev)/max(1,norm(xPrev)) <= TolX 
      Flag    = FLAG_TOLX;
      Message = MESSAGE_TOLX;
      break
    elseif f <= min(fPrev) && abs(min(fPrev)-f)/max(1,abs(fPrev(end))) <= TolFun
      Flag    = FLAG_TOLFUN;
      Message = MESSAGE_TOLFUN;
      break
      
    % Check function evaluation/iteration cap
    elseif iter >= MaxIter 
      Flag    = FLAG_MAXITER;
      Message = MESSAGE_MAXITER;
      break
    elseif FunEvals >= MaxFunEvals
      Flag    = FLAG_MAXFUNEVALS;
      Message = MESSAGE_MAXFUNEVALS;
      break
    end
  end
  
  % -------------------- Cleanup and exit --------------------
  
  Trace.f        = Trace.f(1:iter+1);
  Trace.FunEvals = Trace.FunEvals(1:iter+1);
  if CheckOpt
    Trace.Optimality = Trace.Optimality(1:iter+1);
  end
  
  if Display && mod(iter,PrintEvery) > 0
    if CheckOpt
      fprintf(' %4d  %8d  %12.4e  %12.4e  %12.4e\n',...
        iter, FunEvals, step, f, opt);
    else
      fprintf(' %4d  %8d  %12.4e  %12.4e\n',...
        iter, FunEvals, step, f);
    end
  end
  
  output = struct(...
    'Flag'      , Flag            ,...
    'FunEvals'  , FunEvals        ,...
    'Iterations', iter            ,...
    'Trace'     , Trace            ...
    );
  if CheckOpt
     output.Optimality = opt;
  end
  
  if Display
    if CheckOpt
      fprintf(' %s\n',repmat('-',1,57));
      fprintf(' %s\n',Message)
      fprintf(' %s\n',repmat('-',1,57));
    else
      fprintf(' %s\n',repmat('-',1,43));
      fprintf(' %s\n',Message)
      fprintf(' %s\n',repmat('-',1,43));
    end
  end
