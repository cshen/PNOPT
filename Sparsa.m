function [x, f, output] = Sparsa(smoothF, nonsmoothF, x, options)
% Sparsa : Structured reconstruction by separable approximation
% 
% [x, f, output] = Sparsa(smoothF, nonsmoothF, x) starts at x and seeks a minimizer
%   of the objective function in composite form. smoothF is a handle to a function
%   that returns the smooth function value and gradient. nonsmoothF is a handle 
%   to a function that returns the nonsmooth function value and proximal mapping.   
% 
% [x, f, output] = Sparsa(smoothF, nonsmoothF, x, options) replaces the default 
%   optimization options with those in options, a structure created using the 
%   PNoptimset function.
% 
  REVISION = '$Revision: 0.4.2$';
  DATE     = '$Date: July 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Parse options ============
  
  SparsaOptions = PNoptimset(...
    'checkOpt'      , 1      ,... % Check opt (requires prox evaluation)
    'debug'         , 0      ,... % debug mode 
    'descParam'     , 0.0001 ,... % sufficient descent parameter
    'display'       , 100    ,... % display frequency (<= 0 for no display) 
    'lineSearchMem' , 10     ,... % number of previous function values to save
    'maxfunEvals'   , 50000  ,... % max number of function evaluations
    'maxIter'       , 5000   ,... % max number of iterations
    'funTol'        , 1e-9   ,... % stopping tolerance on objective function 
    'optTol'        , 1e-6   ,... % stopping tolerance on opt
    'xtol'          , 1e-9    ... % stopping tolerance on solution
    );
  
  if nargin > 3
    options = PNoptimset(SparsaOptions, options);
  else
    options = SparsaOptions;
  end
  
  checkOpt      = options.checkOpt;
  debug         = options.debug;
  descParam     = options.descParam;
  display       = options.display;
  lineSearchMem = options.lineSearchMem;
  maxfunEvals   = options.maxfunEvals;
  maxIter       = options.maxIter;
  funTol        = options.funTol;
  optTol        = options.optTol;
  xtol          = options.xtol;
  
% ============ Initialize variables ============
  
  FLAG_OPT         = 1;
  FLAG_XTOL        = 2;
  FLAG_FUNTOL      = 3;
  FLAG_MAXITER     = 4;
  FLAG_MAXFUNEVALS = 5;
  
  MSG_OPT         = 'Optimality below optTol.';
  MSG_XTOL        = 'Relative change in x below xtol.';
  MSG_FUNTOL      = 'Relative change in function value below funTol.';
  MSG_MAXITER     = 'Max number of iterations reached.';
  MSG_MAXFUNEVALS = 'Max number of function evaluations reached.';
  
  iter = 0; 
  loop = 1;
  
  Trace.f         = zeros(maxIter+1,1);
  Trace.funEvals  = zeros(maxIter+1,1);
  Trace.proxEvals = zeros(maxIter+1,1);
  if checkOpt
    Trace.opt     = zeros(maxIter+1,1);
  end
  
  if debug
    Trace.normSearchDir   = zeros(maxIter,1);
    Trace.lineSearchFlag  = zeros(maxIter,1);
    Trace.lineSearchIters = zeros(maxIter,1);
  end
  
  if display > 0    
    if checkOpt
      fprintf(' %s\n',repmat('=',1,64));
      fprintf('                   SPG v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,64));
      fprintf(' %4s   %6s  %6s  %12s  %12s  %12s \n',...
        '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'Optimality');
      fprintf(' %s\n',repmat('-',1,64));
    else
      fprintf(' %s\n',repmat('=',1,50));
      fprintf('            SPG v.%s (%s)\n', REVISION, DATE);
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
    [~, xProx ] = nonsmoothF( x - Df ,1);
         opt    = norm( xProx - x ,'inf');
  end
  
  Trace.f(1)         = f;
  Trace.funEvals(1)  = funEvals;
  Trace.proxEvals(1) = proxEvals;
  if checkOpt
    Trace.opt(1)     = opt; 
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
  
  if checkOpt && opt <= optTol
    flag    = FLAG_OPT;
    message = MSG_OPT;
    loop    = 0;
  end

% ============ Main Loop ============
  
  while loop
    iter = iter+1; 
    
  % ------------ Compute search direction ------------
    
    if iter > 1
      s  = x-xPrev;
      y  = Df-DfPrev;
      BBstep  = (y'*s)/(y'*y);
      if BBstep <= 1e-9 || 1e9 <= BBstep
        BBstep = min(1,1/norm(Df,1));
      end
    else
      BBstep = min(1,1/norm(Df,1));
    end
    
  % ------------ Conduct line search ------------
    
    xPrev   = x;
    if iter+1 > lineSearchMem
      fPrev = [fPrev(2:end), f];
    else
      fPrev(iter) = f;
    end
    DfPrev  = Df;
    
    [x, f, Df, step, lineSearchFlag ,lineSearchIters] = ...
      CurvySearch(x, -Df, BBstep, fPrev, -norm(Df)^2, smoothF, nonsmoothF,...
        descParam, xtol, maxfunEvals - funEvals); 
    
  % ------------ Collect data and display status ------------
    
    funEvals  = funEvals + lineSearchIters;
    proxEvals = proxEvals + lineSearchIters;
    if checkOpt
      [~, xProx ] = nonsmoothF( x - Df ,1);
           opt    = norm( xProx - x ,'inf');
    end
    
    Trace.f(iter+1)         = f;
    Trace.funEvals(iter+1)  = funEvals;
    Trace.proxEvals(iter+1) = proxEvals;
    if checkOpt
      Trace.opt(iter+1)     = opt; 
    end
    
    if debug
      Trace.normSearchDir(iter)   = norm(Df);
      Trace.lineSearchFlag(iter)  = lineSearchFlag;
      Trace.lineSearchIters(iter) = lineSearchIters;
    end
    
    if display > 0 && mod(iter,display) == 0
      if checkOpt
        fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f, opt);
      else
        fprintf(' %4d | %6d  %6d  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f);
      end
    end
    
  % ------------ Check stopping criteria ------------
    
    if checkOpt && opt <= optTol
      flag    = FLAG_OPT;
      message = MSG_OPT;
      loop    = 0;
    elseif norm(x-xPrev,'inf')/max(1,norm(xPrev,'inf')) <= xtol 
      flag    = FLAG_XTOL;
      message = MSG_XTOL;
      loop    = 0;
    elseif f <= min(fPrev) && abs(min(fPrev)-f)/max(1,abs(fPrev(end))) <= funTol
      flag    = FLAG_FUNTOL;
      message = MSG_FUNTOL;
      loop    = 0;
    elseif iter >= maxIter 
      flag    = FLAG_MAXITER;
      message = MSG_MAXITER;
      loop    = 0;
    elseif funEvals >= maxfunEvals
      flag    = FLAG_MAXFUNEVALS;
      message = MSG_MAXFUNEVALS;
      loop    = 0;
    end
  end
  
% ============ Cleanup and exit ============
  
  Trace.f         = Trace.f(1:iter+1);
  Trace.funEvals  = Trace.funEvals(1:iter+1);
  Trace.proxEvals = Trace.proxEvals(1:iter+1);
  if checkOpt
    Trace.opt     = Trace.opt  (1:iter+1);
  end
  
  if debug
    Trace.normSearchDir   = Trace.normSearchDir(1:iter);
    Trace.lineSearchFlag  = Trace.lineSearchFlag(1:iter);
    Trace.lineSearchIters = Trace.lineSearchIters(1:iter);
  end
  
  if display > 0 && mod(iter,display) > 0
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
    'iters'     , iter      ,...
    'options'   , options   ,...
    'proxEvals' , proxEvals ,...
    'Trace'     , Trace      ...
    );
  if checkOpt
     output.opt = opt;
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
