function [x, f, output] = QuasiNewton(fun, x, options)
% QuasiNewton : Quasi-Newton Methods
% 
% [x, f, output] = QuasiNewton(fun, x) starts at x and seeks a minimizer of the
%   objective function. fun is a handle to a function that returns the objective
%   function value and gradient.
% 
% [x, f, output] = QuasiNewton(fun, x, options) replaces the default options 
%   with those in options, a struct created using the PNoptimset function.
% 
  REVISION = '$Revision: 0.4.2$';
  DATE     = '$Date: July 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Parse options ============
  
  QNoptions = PNoptimset(...
    'curvParam'   , 0.9     ,... % curvature condition parameter
    'debug'       , 0       ,... % debug mode 
    'descParam'   , 0.0001  ,... % Sufficient descent parameter
    'display'     , 10      ,... % display frequency (<= 0 for no display) 
    'LbfgsMem'    , 50      ,... % Number of L-BFGS corrections
    'maxfunEvals' , 5000    ,... % Max number of function evaluations
    'maxIter'     , 500     ,... % Max number of iterations
    'method'      , 'Bfgs' ,... % method for choosing search directions
    'funTol'      , 1e-9    ,... % Stopping tolerance on objective function 
    'optTol'      , 1e-6    ,... % Stopping tolerance on opt
    'xtol'        , 1e-9     ... % Stopping tolerance on solution
    );
  
  if nargin > 2
    options = PNoptimset(QNoptions, options);
  else
    options = QNoptions;
  end
  
  curvParam    = options.curvParam;
  debug        = options.debug;
  descParam    = options.descParam;
  display      = options.display;
  maxfunEvals  = options.maxfunEvals;
  maxIter      = options.maxIter;
  method       = options.method;
  switch method
    case 'Bfgs'
      evalHess = 0;
    case 'Lbfgs'
      evalHess = 0;
      LbfgsMem = options.LbfgsMem;
    case 'Newton'
      evalHess = 1;
  end
  funTol       = options.funTol;
  optTol       = options.optTol;
  xtol         = options.xtol;
  
  % ============ Initialize local variables ============
  
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
  
  Trace.f        = zeros(maxIter+1,1);
  Trace.funEvals = zeros(maxIter+1,1);
  Trace.opt      = zeros(maxIter+1,1);
  Trace.step     = zeros(maxIter,1);
  
  if debug
    Trace.normSearchDir   = zeros(maxIter,1);
    Trace.gtd             = zeros(maxIter,1);
    Trace.lineSearchFlag  = zeros(maxIter,1);
    Trace.lineSearchIters = zeros(maxIter,1);
  end
  
  if display > 0
    fprintf(' %s\n',repmat('=',1,56));
    fprintf('           QuasiNewton v.%s (%s)\n', REVISION, DATE);
    fprintf(' %s\n',repmat('=',1,56));
    fprintf(' %4s   %6s  %12s  %12s  %12s \n',...
      '','Fun.', 'Step len.', 'Obj. val.', 'Optimality');
    fprintf(' %s\n',repmat('-',1,56));
  end
  
  % ============ Evaluate objective function at starting x ============ 
  
  if evalHess
    [f, Df, Hf] = fun(x);
  else
    [f, Df] = fun(x);
  end
  
  % ============ Start collecting data for display and output ============ 
  
  funEvals = 1;
  opt      = norm( Df ,'inf');
  
  Trace.f(1)        = f;
  Trace.funEvals(1) = funEvals;
  Trace.opt(1)      = opt;
  
  if display > 0
    fprintf(' %4d | %6d  %12s  %12.4e  %12.4e\n',...
      iter, funEvals, '', f, opt);
  end
  
  % ============ Check if starting x is optimal ============ 
  
  if opt <= optTol      
    flag    = FLAG_OPT;
    message = MSG_OPT;
    loop    = 0;
  end
  
% ============ Main Loop ============
  
  while loop
    iter = iter+1; 
    
    % ------------ Compute search direction ------------
    
    switch method
      
      % BFGS method
      case 'Bfgs'
        if iter > 1
          s =  x - xPrev;
          y = Df - DfPrev;
          qty1 = cholB'*(cholB*s);
          if s'*y > 1e-9
            cholB = cholupdate(cholupdate(cholB, y/sqrt(y'*s)), qty1/sqrt(s'*qty1),'-');
          end
          searchDir = -cholB\(cholB'\Df);
        else
          cholB     = eye(length(x));
          searchDir = -Df;
        end
        
      % Limited-memory BFGS method
      case 'Lbfgs'
        if iter > 1
          s =  x - xPrev;
          y = Df - DfPrev;
          if y'*s > 1e-9
            if size(sPrev,2) > LbfgsMem
              sPrev = [sPrev(:,2:LbfgsMem), s];
              yPrev = [yPrev(:,2:LbfgsMem), y];
              de    = (y'*y)/(y'*s);
            else
              sPrev = [sPrev, s]; %#ok<AGROW>
              yPrev = [yPrev, y]; %#ok<AGROW>
              de    = (y'*y)/(y'*s);
            end
          end
          searchDir = LbfgsSearchDir(sPrev, yPrev, de, Df);
        else
          sPrev     = zeros(length(x), 0);
          yPrev     = zeros(length(x), 0);
          de        = 1;
          searchDir = -Df;
        end
        
      % Newton's method
      case 'Newton'
        if isa(Hf,'function_handle')
          searchDir = pcg(Hf, -Df, min(0.5,sqrt(opt))*opt);
        elseif isnumeric(Hf)
          searchDir = Hf\(-Df);
        end
    end
    
    % ------------ Conduct line search ------------
    
    xPrev  = x;
    fPrev  = f;
    DfPrev = Df;
    
    if evalHess
      [x, f, Df, Hf, step, lineSearchFlag, lineSearchIters] = ...
        MoreThuenteSearch(fun, x, f, Df, searchDir, 1, descParam, curvParam,...
          max(xtol,1e-9), maxfunEvals - funEvals);
    else
      if iter > 1
        [x, f, Df, step, lineSearchFlag, lineSearchIters] = ...
          MoreThuenteSearch(fun, x, f, Df, searchDir, 1, descParam, curvParam,...
          max(xtol,1e-9), maxfunEvals - funEvals);
      else
        [x, f, Df, step, lineSearchFlag, lineSearchIters] = ...
          MoreThuenteSearch(fun, x, f, Df, searchDir, min(1,1/norm(Df)),...
          descParam, curvParam, max(xtol,1e-9), maxfunEvals - funEvals);
      end
    end
    
    % ------------ Collect data for display and output ------------
    
    funEvals = funEvals + lineSearchIters;   
    opt      = norm( Df ,'inf');
    
    Trace.f(iter+1)        = f;
    Trace.funEvals(iter+1) = funEvals;
    Trace.opt  (iter+1)    = opt;
    Trace.step(iter)       = step;
    
    if debug
      Trace.normSearchDir(iter)   = norm(searchDir);
      Trace.gtd(iter)             = Df'*searchDir;
      Trace.lineSearchFlag(iter)  = lineSearchFlag;
      Trace.lineSearchIters(iter) = lineSearchIters;
    end
    
    if display > 0 && mod(iter,display) == 0
      fprintf(' %4d | %6d  %12.4e  %12.4e  %12.4e\n',...
        iter, funEvals, step, f, opt);
    end
    
    % ------------ Check stopping criteria ------------
    
    if opt <= optTol
      flag    = FLAG_OPT;
      message = MSG_OPT;
      loop    = 0;
    elseif norm(x-xPrev,'inf')/max(1,norm(xPrev,'inf')) <= xtol 
      flag    = FLAG_XTOL;
      message = MSG_XTOL;
      loop    = 0;
    elseif f <= fPrev && (fPrev-f)/max(1,abs(fPrev)) <= funTol
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
  
  Trace.f        = Trace.f(1:iter+1);
  Trace.funEvals = Trace.funEvals(1:iter+1);
  Trace.opt      = Trace.opt  (1:iter+1);
  
  if debug
    Trace.normSearchDir   = Trace.normSearchDir(1:iter);
    Trace.gtd             = Trace.gtd(1:iter);
    Trace.lineSearchFlag  = Trace.lineSearchFlag(1:iter);
    Trace.lineSearchIters = Trace.lineSearchIters(1:iter);
  end
  
  if display > 0 && mod(iter,display) > 0
    fprintf(' %4d | %6d  %12.4e  %12.4e  %12.4e\n',...
      iter, funEvals, step, f, opt);
  end
      
  output = struct(...
    'flag'     , flag     ,...
    'funEvals' , funEvals ,...
    'iters'    , iter     ,...
    'opt'      , opt    ,...
    'options'  , options  ,...
    'Trace'    , Trace     ...
    );
  
  if display > 0
    fprintf(' %s\n',repmat('-',1,56));
    fprintf(' %s\n',message)
    fprintf(' %s\n',repmat('-',1,56));
  end
