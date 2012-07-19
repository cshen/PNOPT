function [x, f, output] = ProxNewtonDual(smoothF, nonsmoothF, conjF, x, options)
% ProxNewtonDual : Dual Proximal Newton-type methods
% 
% [x, f, output] = ProxNewtonDual(smoothF, nonsmoothF, conjF, x) starts at x and 
%   seeks a minimizer of the objective function in composite form. smoothF is a 
%   handle to a function that returns the smooth function value and gradient. 
%   nonsmoothF and conjF are handles to functions that returns the value and 
%   proximal mapping of the nonsmooth function and its conjugate.
% 
% [x, f, output] = ProxNewtonDual(smoothF, nonsmoothF, x, options) replaces the   
%   default optimization options with those in options, a structure created using
%  the PNoptimset function.
% 
  REVISION = '$Revision: 0.1.0$';
  DATE     = '$Date: July 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Process options ============
  
  defaultOptions = PNoptimset(...
    'debug'       , 0             ,... % debug mode 
    'descParam'   , 0.0001        ,... % Armijo condition parameter
    'display'     , 10            ,... % display frequency (<= 0 for no display) 
    'LbfgsMem'    , 50            ,... % number of L-BFGS corrections
    'maxfunEvals' , 5000          ,... % max number of function evaluations
    'maxIter'     , 500           ,... % max number of iterations
    'method'      , 'Lbfgs'       ,... % method for choosing search directions
    'subMethod'   , 'QuasiNewton' ,... % solver for solving subproblems
    'funTol'      , 1e-9          ,... % stopping tolerance on relative change in the objective function 
    'optTol'      , 1e-6          ,... % stopping tolerance on opt
    'xtol'        , 1e-9           ... % stopping tolerance on solution
    );
  
  minFunc_options.display = 0;
  
  if nargin > 4
    options = PNoptimset(defaultOptions, options);
  else
    options = defaultOptions;
  end
  
  debug       = options.debug;
  descParam   = options.descParam;
  display     = options.display;
  maxfunEvals = options.maxfunEvals;
  maxIter     = options.maxIter;
  method      = options.method;
  switch method
    case 'Bfgs'
      evalHess = 0;
    case 'Lbfgs'
      evalHess = 0;
      LbfgsMem = options.LbfgsMem;
    case 'Newton'
      evalHess = 1;
  end
  subMethod   = options.subMethod;
  funTol      = options.funTol;
  optTol      = options.optTol;
  xtol        = options.xtol;
  
% ------------ Set subproblem solver options ------------
  
  switch subMethod
    case 'minfunc'
      if isfield(options, 'minFunc_options') && ~isempty(options.minFunc_options)
        minFunc_options = mergestruct(minFunc_options, options.minFunc_options);
      end
  end
  
% ============ Initialize variables ============
  
  FLAG_OPT         = 1;
  FLAG_XTOL        = 2;
  FLAG_FUNTOL      = 3;
  FLAG_MAXITER     = 4;
  FLAG_MAXFUNEVALS = 5;
  FLAG_OTHER       = 6;
  
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
  Trace.opt       = zeros(maxIter+1,1);
  
  if debug
    Trace.dualityGap      = zeros(maxIter,1);
    Trace.normSearchDir   = zeros(maxIter,1);
    Trace.lineSearchFlag  = zeros(maxIter,1);
    Trace.lineSearchIters = zeros(maxIter,1);
    Trace.subFlags = zeros(maxIter,1);
    Trace.subIters = zeros(maxIter,1);
    Trace.subopt   = zeros(maxIter,1);
  end
  
  if display > 0  
    fprintf(' %s\n',repmat('=',1,64));
    fprintf('             ProxNewtonDual  v.%s (%s)\n', REVISION, DATE);
    fprintf(' %s\n',repmat('=',1,64));
    fprintf(' %4s   %6s  %6s  %12s  %12s  %12s \n',...
      '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'Optimality');
    fprintf(' %s\n',repmat('-',1,64));
  end
  
% ------------ Evaluate objective function at starting x ------------
  
  if evalHess
    [f, Df, Hf] = smoothF(x);
  else
    [f, Df] = smoothF(x);
  end
   h = nonsmoothF(x);
   f = f + h;
  
% ------------ Start collecting data for display and output ------------
  
  funEvals  = 1;
  proxEvals = 0;
  [~, xProx ] = nonsmoothF( x - Df, 1);
       opt    = norm( xProx - x ,'inf');
  
  Trace.f(1)         = f;
  Trace.funEvals(1)  = funEvals;
  Trace.proxEvals(1) = proxEvals;
  Trace.opt(1)       = opt; 
  
  if display > 0
    fprintf(' %4d | %6d  %6d  %12s  %12.4e  %12.4e\n',...
      iter, funEvals, proxEvals, '', f, opt);
  end
  
% ------------ Check if starting x is optimal ------------
  
  if opt <= optTol
    flag    = FLAG_OPT;
    message = MSG_OPT;
    loop    = 0;
  end

% ============ Main Loop ============
  
  while loop
    iter = iter+1; 
    
  % ------------ Update Hessian approximation ------------
    
    switch method
      
      % BFGS method
      case 'Bfgs'
        if iter > 1
          s =  x - xPrev;
          y = Df - DfPrev;
          qty1 = R'*(R*s);
          if s'*y > 1e-9
            R = cholupdate(cholupdate(R, y/sqrt(y'*s)), qty1/sqrt(s'*qty1),'-');
          end
          Hf      = @(x) R'*(R*x);
          Bf      = @(x) R\(R'\x);
          proxCtr = x - R\(R'\Df);
        else
          R  = eye(length(x));
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
              et    = (y'*y)/(y'*s);
            else
              sPrev = [sPrev, s]; %#ok<AGROW>
              yPrev = [yPrev, y]; %#ok<AGROW>
              et    = (y'*y)/(y'*s);
            end
          end
          Hf      = LbfgsProd(sPrev, yPrev, et);
          Bf      = @(x) LbfgsSearchDir(sPrev, yPrev, et, -x);
          proxCtr = x - LbfgsSearchDir(sPrev, yPrev, et, Df);
        else
          sPrev = zeros(length(x), 0);
          yPrev = zeros(length(x), 0);
        end
        
      % Newton's method
      case 'Newton'
        if iter > 1
          s  =  x - xPrev;
          y  = Df - DfPrev;
          et = (y'*y)/(y'*s);
        end
        
        if isa(Hf,'function_handle')
          newtonDir = pcg(Hf, -Df, min(0.5,sqrt(opt))*opt);
        elseif isnumeric(Hf)
          newtonDir = Hf\(-Df);
        end
        proxCtr = x + newtonDir;
    end
    
  % ------------ Solve subproblem for a search direction ------------
    
    if iter > 1 || evalHess      
      subDualF = @(v) pnopt_subDual(mu, Bf, proxCtr, v, nonsmoothF, conjF);
      
      switch subMethod
        
        % minFunc
        case 'minfunc'          
          mu = 1/et;
                    
          [ xProxDual, subDualf, minFunc_flag, minFunc_output ] = ...
            minFunc(subDualF, proxCtr, minFunc_options);  
          
        % ------------ Collect data from subproblem solve ------------
          
          subIters     = minFunc_output.iterations;
          subProxEvals = minFunc_output.funcCount;
          
          if debug
            switch minFunc_flag
              case 0
                switch minFunc_output.message
                  case 'Reached Maximum Number of Function Evaluations'
                    subFlag = FLAG_MAXFUNEVALS;
                  case 'Reached Maximum Number of Iterations'
                    subFlag = FLAG_MAXITER;
                end
              case 1
                subFlag     = FLAG_OPT;
              case 2
                switch minFunc_output.message
                  case 'Step Size below progTol'
                    subFlag = FLAG_XTOL;
                  case 'Function Value changing by less than progTol'
                    subFlag = FLAG_FUNTOL;
                end 
              otherwise
                subFlag = FLAG_OTHER;
            end
            subopt = minFunc_output.firstorderopt;
          end
      end
      [ ~, xProx ] = nonsmoothF( proxCtr - xProxDual/mu, 1/mu ); 
      searchDir = xProx -x;
    else
      subIters     = 0;
      subProxEvals = 0;
      
      if debug 
        subFlag  = 0;
        subopt = 0;
      end
      
      searchDir = -Df;
    end
    
  % ------------ Conduct line search ------------
    
    xPrev  = x;
    fPrev  = f;
    DfPrev = Df;
    
    if evalHess
      [x, f, Df, Hf, step, lineSearchFlag ,lineSearchIters] = ...
        LineSearch(x, searchDir, 1, f, h, Df'*searchDir, smoothF, nonsmoothF,...
          descParam, xtol, maxfunEvals - funEvals); 
    else
      if iter > 1
        [x, f, Df, step, lineSearchFlag ,lineSearchIters] = ...
          LineSearch(x, searchDir, 1, f, h, Df'*searchDir, smoothF, nonsmoothF,...
            descParam, xtol, maxfunEvals - funEvals);
      else
        [x, f, Df, step, lineSearchFlag ,lineSearchIters] = ...
          CurvySearch(x, searchDir, max(min(1,1/norm(Df)), xtol), f, Df'*searchDir,...
            smoothF, nonsmoothF, descParam, xtol, maxfunEvals - funEvals); 
      end
    end 
    
  % ------------ Collect data for display and output ------------
    
    funEvals  =  funEvals + lineSearchIters;
    proxEvals = proxEvals + lineSearchIters + subProxEvals;
    [~, xProx ] = nonsmoothF( x - Df, 1);
         opt    = norm( xProx - x ,'inf');
    
    Trace.f(iter+1)         = f;
    Trace.funEvals(iter+1)  = funEvals;
    Trace.proxEvals(iter+1) = proxEvals;
    Trace.opt(iter+1)       = opt;
    
    if debug
      Trace.normSearchDir(iter)   = norm(searchDir);
      Trace.lineSearchFlag(iter)  = lineSearchFlag;
      Trace.lineSearchIters(iter) = lineSearchIters;
      Trace.subFlags(iter) = subFlag;
      Trace.subIters(iter) = subIters;
      Trace.subopt(iter)   = subopt;
    end
    
    if display > 0 && mod(iter,display) == 0
      fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f, opt);
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
    elseif abs(fPrev-f)/max(1,abs(fPrev)) <= funTol
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
  Trace.opt       = Trace.opt  (1:iter+1);
  
  if debug
    Trace.normSearchDir   = Trace.normSearchDir(1:iter);
    Trace.lineSearchFlag  = Trace.lineSearchFlag(1:iter);
    Trace.lineSearchIters = Trace.lineSearchIters(1:iter);
    Trace.subFlags = Trace.subFlags(1:iter);
    Trace.subIters = Trace.subIters(1:iter);
    Trace.subopt   = Trace.subopt(1:iter);
  end
  
  if display > 0 && mod(iter,display) > 0
    fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
      iter, funEvals, proxEvals, step, f, opt);
  end
  
  output = struct(...
    'flag'      , flag      ,...
    'funEvals'  , funEvals  ,...
    'iters'     , iter      ,...
    'opt'       , opt       ,...
    'options'   , options   ,...
    'proxEvals' , proxEvals ,...
    'Trace'     , Trace      ...
    );
  
  if display > 0
    fprintf(' %s\n',repmat('-',1,64));
    fprintf(' %s\n',message)
    fprintf(' %s\n',repmat('-',1,64));
  end
  
  clear global TfocsDf TfocsOptim
  