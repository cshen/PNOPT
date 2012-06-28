function [x, f, output] = ProxNewton(smoothF, nonsmoothF, x, varargin)
% ProxNewton : Proximal Newton-type methods
% 
% [x, f, output] = ProxNewton(smoothF, nonsmoothF, x) starts at x and seeks a
%   minimizer of the objective function. smoothF is a handle to a function that
%   returns the smooth function value and gradient. nonsmoothF is a handle to a
%   function that returns the nonsmooth function value and proximal point.
% 
% [x, f, output] = ProxNewton(smoothF, nonsmoothF, x, options) replaces the   
%   default optimization options replaced with values in options, a structure
%   created using the SetPNoptOptions function.
% 
  REVISION = '$Revision: 0.2.0$';
  DATE     = '$Date: June 24, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ------------ Initialize ------------
  
  % Set default options for subproblem
  spgOptions = SetPNoptOptions(...
    'checkOpt'   , 1 ,... 
    'display'    , 0 ,...
    'maxfunEvals', 500 ,...
    'maxIter'    , 50  ,...
    'TolOpt'     , 1e-3 ...
    );
  
  TfocsOpts = struct(...
    'alg'       , 'N83',...
    'maxIts'    , 50 ,...
    'printEvery', 0  ,...
    'restart'   , -Inf ,...
    'tol'       , 1e-3  ...
    );

  QNoptions = SetPNoptOptions(...
    'display'    , 1    ,...
    'maxfunEvals', 5000 ,...
    'maxIter'    , 500 ,...
    'TolFun'     , 0    ...
    );
  
  % Set default options
  defaultOptions = SetPNoptOptions(...
    'checkOpt'         , 1         ,... % Check optimality (requires prox evaluation)
    'display'          , 1         ,... % display level 
    'LbfgsCorrections' , 20        ,... % Number of L-BFGS corrections
    'maxfunEvals'      , 5000      ,... % Max number of function evaluations
    'maxIter'          , 500       ,... % Max number of iterations
    'method'           , 'Lbfgs'   ,... % method for choosing search directions
    'printEvery'       , 10        ,... % display frequency (iterations)
    'spgOptions'       , spgOptions,... % Options for solving subproblems using spg
    'subproblemMethod' , 'Tfocs'   ,... % Solver for solving subproblems
    'TfocsOpts'        , TfocsOpts ,... % Options for solving subproblems using TFOCS
    'TolFun'           , 1e-9      ,... % Stopping tolerance on objective function 
    'TolOpt'           , 1e-6      ,... % Stopping tolerance on optimality
    'TolX'             , 1e-9       ... % Stopping tolerance on solution
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
    options = SetPNoptOptions(defaultOptions, varargin{1});
  else
    options = defaultOptions;
  end
  
  checkOpt          = options.checkOpt;
  display           = options.display;
  LbfgsCorrections  = options.LbfgsCorrections;
  maxfunEvals       = options.maxfunEvals;
  maxIter           = options.maxIter;
  method            = options.method;
  printEvery        = options.printEvery;
  subproblemMethod  = options.subproblemMethod;
  TolFun            = options.TolFun;
  TolOpt            = options.TolOpt;
  TolX              = options.TolX;
  
  if strcmp(method,'Newton')
    Hess = options.Hess;
  end
  
  switch subproblemMethod
    case 'spg'
      spgOptions    = options.BbOptions;
      SubproblemTol = spgOptions.TolOpt;
    case 'Tfocs'
      TfocsOpts     = options.TfocsOpts;
      SubproblemTol = TfocsOpts.tol;
  end
  
  iter            = 0; 
  Trace.f         = zeros(maxIter+1,1);
  Trace.funEvals  = zeros(maxIter+1,1);
  Trace.proxEvals = zeros(maxIter+1,1);
  if checkOpt
    Trace.optimality = zeros(maxIter+1,1);
  end
  
  if display    
    if checkOpt
      fprintf(' %s\n',repmat('=',1,64));
      fprintf('                ProxNewton v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,64));
      fprintf(' %4s   %6s  %6s  %12s  %12s  %12s \n',...
        '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'optimality');
      fprintf(' %s\n',repmat('-',1,64));
    else
      fprintf(' %s\n',repmat('=',1,50));
      fprintf('         ProxNewton v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,50));
      fprintf(' %4s   %6s  %6s  %12s  %12s \n',...
        '','Fun.', 'Prox', 'Step len.', 'Obj. val.');
      fprintf(' %s\n',repmat('-',1,50));
    end
  end
  
  % ------------ Evaluate objective function at starting x ------------ 
  
  [f, Df] = smoothF(x);
   h      = nonsmoothF(x);
   f      = f + h;
  if strcmp(method,'Newton')  % Evaluate Hessian if necessary.
    Hf    = Hess(x);
    if ~isa(Hf,'numeric') && ~isa(Hf,'function_handle')  
      error('ProxNewton:BadHess', 'options.Hess must be a function handle or numeric array')
    end
  end
  
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
      fprintf(' %4d | %6d  %6d  %12s  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, '', f, opt);
    else
      fprintf(' %4d | %6d  %6d  %12s  %12.4e \n',...
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
          Hf = @(x) R'*(R*x);
              
        else
          R  = eye(length(x));
        end
        
        if strcmp(subproblemMethod,'smoothed')
          Be = @(x) R\(R'\x);
          Dx = -R\(R'\Df);
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
          Hf = LbfgsProd(sPrev, yPrev, et);
        else
          sPrev = zeros(length(x), 0);
          yPrev = zeros(length(x), 0);
          et    = 1;
        end
        
        if strcmp(subproblemMethod,'smoothed')
          Be = @(x) LbfgsSearchDir(sPrev, yPrev, et, -x);
          Dx = LbfgsSearchDir(sPrev, yPrev, et, Df);
        end
        
      % Newton's method
      case 'Newton'        
        if strcmp(subproblemMethod,'smoothed')
          Be = @(x) Hf\x;
          Dx = -Hf\Df;
        end
    end
    
    % Solve subproblem to obtain search direction
    if iter > 1 || strcmp(method,'Newton')
      switch subproblemMethod
        case 'spg'
          spgOptions = SetPNoptOptions(spgOptions,...
            'TolOpt', SubproblemTol ...
            );
        
          [z, ~, spgOutput] = ...
            spg(@(z) QuadF(Hf, Df, x ,z), nonsmoothF, x, spgOptions);
          
          proxEvals = proxEvals + spgOutput.proxEvals;
          
          % If SPG stops early, then make stopping tolerance smaller.     
          if spgOutput.Iterations < spgOptions.maxIter  
            SubproblemTol  = max(0.5*SubproblemTol, TolOpt);     
          end
          
        case 'smoothed' % doesn't fucking work
          
        case 'Tfocs'
          TfocsOpts.tol = SubproblemTol;
          
          [z, TfocsOut] = ...
            tfocs(@(z) QuadF(Hf, Df, x ,z), [], nonsmoothF, x, TfocsOpts);
          
          if isfield(TfocsOpts, 'countOpts') && TfocsOpts.countOpts
            proxEvals = proxEvals + TfocsOut.counts(end,6);
          else
            proxEvals = proxEvals + TfocsOut.niter;
          end
          
          if TfocsOut.niter < TfocsOpts.maxIts  
            SubproblemTol = max(0.5*SubproblemTol, TolX);
          end
      end
    end
    
    % ------------ Conduct line search ------------    
    xPrev  = x;
    fPrev  = f;
    DfPrev = Df;
    
    % Conduct line search for a step length that safisfies the Armijo condition
    if iter > 1 
      [x, f, h, Df, step, LSflag ,LSiter] = ...
        LineSearch(x, z-x, 1, f, h, Df'*(z-x), smoothF, nonsmoothF,...
        TolX, maxfunEvals - funEvals); %#ok<ASGLU>
    else
      [x, f, Df, step, LSflag ,LSiter] = ...
        CurvySearch(x, -Df, min(1,1/norm(Df,1)), f, -norm(Df)^2, smoothF, nonsmoothF,...
        TolX, maxfunEvals - funEvals); %#ok<ASGLU>          
    end 
    
    if strcmp(method,'Newton')  % Evaluate Hessian if necessary.
      Hf = Hess(x);
      if ~isa(Hf,'numeric') && ~isa(Hf,'function_handle') 
        error('ProxNewton:BadHess', 'options.Hess must be a function handle or numeric array')
      end
    end
    
    % ------------ Collect data for display and output ------------
    
    funEvals  = funEvals + LSiter;
    proxEvals = proxEvals + LSiter;
    if checkOpt
      [h1, x1]  = nonsmoothF(x-Df,1); %#ok<ASGLU>
      proxEvals = proxEvals + 1;
      opt       = norm(x1-x,'inf'); 
    end
    
    Trace.f(iter+1)         = f;
    Trace.funEvals(iter+1)  = funEvals;
    Trace.proxEvals(iter+1) = proxEvals;
    if checkOpt
      Trace.optimality(iter+1) = opt; 
    end
    
    if display && mod(iter,printEvery) == 0
      if checkOpt
        fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
          iter, funEvals, proxEvals, step, f, opt);
      else
        fprintf(' %4d | %6d  %6d  %12.4e  %12.4e\n',...
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
    elseif abs(fPrev-f)/max(1,abs(fPrev)) <= TolFun
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
      fprintf(' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f, opt);
    else
      fprintf(' %4d | %6d  %6d  %12.4e  %12.4e\n',...
        iter, funEvals, proxEvals, step, f);
    end
  end
  
  output = struct(...
    'flag'      , flag     ,...
    'funEvals'  , funEvals ,...
    'Iterations', iter     ,...
    'method'    , method   ,...
    'options'   , options  ,...
    'proxEvals' , proxEvals,...
    'Trace'     , Trace     ...
    );
  if checkOpt
     output.Optimality = opt;
  end
  
  if display
    if checkOpt
      fprintf(' %s\n',repmat('-',1,64));
      fprintf(' %s\n',Message)
      fprintf(' %s\n',repmat('-',1,64));
    else
      fprintf(' %s\n',repmat('-',1,50));
      fprintf(' %s\n',Message)
      fprintf(' %s\n',repmat('-',1,50));
    end
  end
  