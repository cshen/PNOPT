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
  REVISION = '$Revision: 0.1.2$';
  DATE     = '$Date: June 12, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% -------------------- Initialize --------------------
  
  n = length(x);

  % Set default options for subproblem
  SpgOptions = SetPNoptOptions(...
    'CheckOpt'   , 1   ,... 
    'Display'    , 0   ,...
    'MaxFunEvals', 1000,...
    'MaxIter'    , 100 ,...
    'TolOpt'     , 1e-3 ...
    );
  
  TfocsOpts = struct(...
    'alg'       , 'N83',...
    'maxIts'    , 100  ,...
    'printEvery', 0    ,...
    'restart'   , -Inf ,...
    'tol'       , 1e-4  ...
    );
  
  % Set default options
  DefaultOptions = SetPNoptOptions(...
    'SpgOptions'       , SpgOptions,... % Options for solving subproblems using spg
    'CheckOpt'         , 1         ,... % Check optimality (requires prox evaluation)
    'Display'          , 1         ,... % Display level 
    'LbfgsCorrections' , 10        ,... % Number of L-BFGS corrections
    'MaxFunEvals'      , 5000      ,... % Max number of function evaluations
    'MaxIter'          , 500       ,... % Max number of iterations
    'Method'           , 'Lbfgs'   ,... % Method for choosing search directions
    'PrintEvery'       , 10        ,... % Display frequency (iterations)
    'SmallFirstStep'   , 1         ,... % Choose small step length during first iteration
    'SubproblemMethod' , 'Tfocs'   ,... % Solver for solving subproblems
    'TfocsOpts'        , TfocsOpts ,... % Options for solving subproblems using TFOCS
    'TolFun'           , 1e-9      ,... % Stopping tolerance on objective function 
    'TolOpt'           , 1e-6      ,... % Stopping tolerance on optimality
    'TolX'             , 1e-9      ,... % Stopping tolerance on solution
    'UseHess'          , 0         ,... % Evaluate Hessian/inexact Hessian
    'UseMex'           , 0          ... % Use MEX functions
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
    options = SetPNoptOptions(DefaultOptions, varargin{1});
  else
    options = DefaultOptions;
  end
  
  CheckOpt          = options.CheckOpt;
  Display           = options.Display;
  LbfgsCorrections  = options.LbfgsCorrections;
  MaxFunEvals       = options.MaxFunEvals;
  MaxIter           = options.MaxIter;
  Method            = options.Method;
  PrintEvery        = options.PrintEvery;
  SmallFirstStep    = options.SmallFirstStep;
  SubproblemMethod  = options.SubproblemMethod;
  TolFun            = options.TolFun;
  TolOpt            = options.TolOpt;
  TolX              = options.TolX;
  UseHess           = options.UseHess;
  UseMex            = options.UseMex;
  
  switch SubproblemMethod
    case 'Spg'
      SpgOptions    = options.BbOptions;
      SubproblemTol = SpgOptions.TolOpt;
    case 'Tfocs'
      TfocsOpts     = options.TfocsOpts;
      SubproblemTol = TfocsOpts.tol;
  end
  
  iter = 0; 
  
  % Evaluate objective function at starting x
  [f, Df] = smoothF(x);
   h      = nonsmoothF(x);
   f      = f + h;
  if UseHess  % Evaluate Hessian if necessary.
    Hess = options.Hess;
    if isa(Hess,'double') || isa(Hess,'function_handle')
      Hf = Hess(x);
    else
       error('ProxNewton:BadHess', 'options.Hess is a function handle or double array.')
    end
  end
  
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
  SubproblemOutputs = cell(MaxIter,1);
  
  Trace.f(1)        = f;
  Trace.FunEvals(1) = FunEvals;
  if CheckOpt
    Trace.Optimality(1) = opt; 
  end
  
  if Display    
    if CheckOpt
      fprintf(' %s\n',repmat('=',1,57));
      fprintf('            ProxNewton v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,57));
      fprintf(' %4s  %8s  %12s  %12s  %12s \n',...
        '','F evals', 'Step len.', 'Obj. val.', 'Optimality');
      fprintf(' %s\n',repmat('-',1,57));
      fprintf(' %4d  %8d  %12s  %12.4e  %12.4e\n',...
        iter, FunEvals, '', f, opt);
    else
      fprintf(' %s\n',repmat('=',1,43));
      fprintf('     ProxNewton v.%s (%s)\n', REVISION, DATE);
      fprintf(' %s\n',repmat('=',1,43));
      fprintf(' %4s  %8s  %12s  %12s \n',...
        '','F evals', 'Step len.', 'Obj. val.');
      fprintf(' %s\n',repmat('-',1,43));
      fprintf(' %4d  %8d  %12s  %12.4e \n',...
        iter, FunEvals, '', f);
    end
  end
  
  % Check if starting x is optimal
  if CheckOpt && opt <= TolOpt
    output = struct(...
      'Flag'      , FLAG_OPTIMAL,...
      'FunEvals'  , FunEvals    ,...
      'Iterations', iter        ,...
      'Method'    , Method      ,...
      'Optimality', opt         ,...
      'Options'   , options     ,...
      'Trace'     , Trace        ...
      );
    return
  end

% ------------ Main Loop -------------
  
  while 1
    iter = iter+1; 
    
    % ------------ Compute search direction ------------
    
    switch Method
      % BFGS method
      case 'Bfgs'
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
        
      % Newton's method
      case 'Newton'        
        if isa(Hf,'double')
          Hf = @(x) Hf*x;
        end
    end
    
    % Solve subproblem to obtain search direction
    if iter > 1 || strcmp(Method,'Newton')
      switch SubproblemMethod
        % SPG
        case 'Spg'
          SpgOptions = SetSpgOptions(SpgOptions,...
            'TolOpt', SubproblemTol ...
            );
          [z, ~, SubproblemOutput] = ...
            Spg(@(z) QuadF(Hf, Df, x ,z), nonsmoothF, x, SpgOptions);
          
          % If BB stops early, then make stopping tolerance smaller.     
          if SubproblemOutput.Iterations < SpgOptions.MaxIter  
            SubproblemTol  = max(0.5*SubproblemTol, TolOpt);     
          end
          
        % Tfocs
        case 'Tfocs'
          TfocsOpts.tol = SubproblemTol;
          [z, SubproblemOutput] = ...
            tfocs(@(z) QuadF(Hf, Df, x ,z), [], nonsmoothF, x, TfocsOpts);
          
          if SubproblemOutput.niter < TfocsOpts.maxIts  
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
      [x, f, h, Df, step, LineSearchFlag ,LineSearchFunEvals] = ...
        LineSearch(x, z-x, 1, f, h, Df'*(z-x), smoothF, nonsmoothF,...
        TolX, MaxFunEvals-FunEvals); %#ok<ASGLU>
    else
      if SmallFirstStep 
        [x, f, Df, step, LineSearchFlag ,LineSearchFunEvals] = ...
          CurvySearch(x, -Df, max(min(1,1/norm(Df,1)),TolX), f, -norm(Df)^2, smoothF, nonsmoothF,...
          TolX, MaxFunEvals-FunEvals); %#ok<ASGLU>
      else
        [x, f, Df, step, LineSearchFlag ,LineSearchFunEvals] = ...
          CurvySearch(x, -Df, 1, f, -norm(Df)^2, smoothF, nonsmoothF,...
          TolX, MaxFunEvals-FunEvals); %#ok<ASGLU>          
      end
      
      SubproblemOutput = struct(...
        'Iterations', LineSearchFunEvals ...
        );
    end 
    
    if UseHess
      Hf = Hess(x);
    end
    
    % ------------ Collect data for display and output ------------
    
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
    SubproblemOutputs{iter}  = SubproblemOutput;
    
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
    elseif abs(fPrev-f)/max(1,abs(fPrev)) <= TolFun
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
  Trace.SubproblemOutputs = { SubproblemOutputs{1:iter} };  %#ok<CCAT1>
  
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
    'Flag'      , Flag           ,...
    'FunEvals'  , FunEvals       ,...
    'Iterations', iter           ,...
    'Method'    , Method         ,...
    'Options'   , options        ,...
    'Trace'     , Trace           ...
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
