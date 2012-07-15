function [x, f, output] = ProxNewton(smoothF, nonsmoothF, x, options)
% ProxNewton : Proximal Newton-type methods
% 
% [x, f, output] = ProxNewton(smoothF, nonsmoothF, x) starts at x and seeks a 
%   minimizer of the objective function in composite form. smoothF is a handle
%   to a function that returns the smooth function value and gradient. nonsmoothF
%   is a handle to a function that returns the nonsmooth function value and 
%   proximal mapping. 
% 
% [x, f, output] = ProxNewton(smoothF, nonsmoothF, x, options) replaces the   
%   default optimization options with those in options, a structure created 
%   using the PNoptimset function.
% 
  REVISION = '$Revision: 0.4.2$';
  DATE     = '$Date: July 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Process options ============
  
  PNoptions = PNoptimset(...
    'debug'       , 0       ,... % debug mode 
    'descParam'   , 0.0001  ,... % Sufficient descent parameter
    'display'     , 10      ,... % display frequency (<= 0 for no display) 
    'LbfgsCorrs'  , 50      ,... % Number of L-BFGS corrections
    'maxfunEvals' , 5000    ,... % Max number of function evaluations
    'maxIter'     , 500     ,... % Max number of iterations
    'method'      , 'Lbfgs' ,... % method for choosing search directions
    'subMethod'   , 'Tfocs' ,... % Solver for solving subproblems
    'funTol'      , 1e-9    ,... % Stopping tolerance on relative change in the objective function 
    'optTol'      , 1e-6    ,... % Stopping tolerance on optimality condition
    'xTol'        , 1e-9     ... % Stopping tolerance on solution
    );
  
  SparsaOptions = PNoptimset(...
    'display'     , 0    ,...
    'maxfunEvals' , 5000 ,...
    'maxIter'     , 500   ...
    );
  
  TfocsOpts = struct(...
    'alg'        , 'N83' ,...
    'maxIts'     , 500   ,...
    'printEvery' , 0     ,...
    'restart'    , -Inf   ...
    );
  
  if nargin > 3
    options = PNoptimset(PNoptions, options);
  else
    options = PNoptions;
  end
  
  debug        = options.debug;
  descParam    = options.descParam;
  display      = options.display;
  LbfgsCorrs   = options.LbfgsCorrs;
  maxfunEvals  = options.maxfunEvals;
  maxIter      = options.maxIter;
  method       = options.method;
  switch method
    case 'Bfgs'
      evalHess = 0;
    case 'Lbfgs'
      evalHess = 0;
    case 'Newton'
      evalHess = 1;
  end
  subMethod    = options.subMethod;
  funTol       = options.funTol;
  optTol       = options.optTol;
  xTol         = options.xTol;
  
  switch subMethod
    case 'Sparsa'
      if isfield(options, 'SparsaOptions') && ~isempty(options.SparsaOptions)
        SparsaOptions = PNoptimset(SparsaOptions, options.SparsaOptions);
      end
    case 'Tfocs'
      if isfield(options, 'TfocsOpts') && ~isempty(options.TfocsOpts)
        TfocsOpts = mergestruct(TfocsOpts, options.TfocsOpts);
      end
           
      if debug
        TfocsOpts.countOps = 1;
        TfocsOpts.errFcn   = @(f,x) tfocs_err();
      end
  end
  
  % ============ Initialize variables ============
  
  FLAG_OPT         = 1;
  FLAG_TOLX        = 2;
  FLAG_TOLFUN      = 3;
  FLAG_MAXITER     = 4;
  FLAG_MAXFUNEVALS = 5;
  FLAG_OTHER       = 6;
  
  MSG_OPT         = 'Optimality below optTol.';
  MSG_TOLX        = 'Relative change in x below xTol.';
  MSG_TOLFUN      = 'Relative change in function value below funTol.';
  MSG_MAXITER     = 'Max number of iterations reached.';
  MSG_MAXFUNEVALS = 'Max number of function evaluations reached.';
  
  iter        = 0; 
  loop        = 1;
  forcingTerm = 0.5;
  
  Trace.f         = zeros(maxIter+1,1);
  Trace.funEvals  = zeros(maxIter+1,1);
  Trace.proxEvals = zeros(maxIter+1,1);
  Trace.opt       = zeros(maxIter+1,1);
  
  if debug
    Trace.normSearchDir   = zeros(maxIter,1);
    Trace.lineSearchFlag  = zeros(maxIter,1);
    Trace.lineSearchIters = zeros(maxIter,1);
    Trace.subFlags = zeros(maxIter,1);
    Trace.subIters = zeros(maxIter,1);
    Trace.subopt   = zeros(maxIter,1);
  end
  
  if display > 0  
    fprintf(' %s\n',repmat('=',1,64));
    fprintf('               ProxNewton  v.%s (%s)\n', REVISION, DATE);
    fprintf(' %s\n',repmat('=',1,64));
    fprintf(' %4s   %6s  %6s  %12s  %12s  %12s \n',...
      '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'Optimality');
    fprintf(' %s\n',repmat('-',1,64));
  end
  
  % ============ Evaluate objective function at starting x ============ 
  
  if evalHess
    [f, Df, Hf] = smoothF(x);
  else
    [f, Df] = smoothF(x);
  end
   h = nonsmoothF(x);
   f = f + h;
  
  % ============ Start collecting data for display and output ============ 
  
  funEvals  = 1;
  proxEvals = 0;
  [~, xProx ] = nonsmoothF( x - Df ,1);
       opt    = norm( xProx - x ,'inf');
  
  Trace.f(1)         = f;
  Trace.funEvals(1)  = funEvals;
  Trace.proxEvals(1) = proxEvals;
  Trace.opt(1)       = opt; 
  
  if display > 0
    fprintf(' %4d | %6d  %6d  %12s  %12.4e  %12.4e\n',...
      iter, funEvals, proxEvals, '', f, opt);
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
          Hf = @(x) R'*(R*x);
        else
          R  = eye(length(x));
        end

      % Limited-memory BFGS method
      case 'Lbfgs'
        if iter > 1
          s =  x - xPrev;
          y = Df - DfPrev;
          if y'*s > 1e-9
            if size(sPrev,2) > LbfgsCorrs
              sPrev = [sPrev(:,2:LbfgsCorrs), s];
              yPrev = [yPrev(:,2:LbfgsCorrs), y];
              de    = (y'*y)/(y'*s);
            else
              sPrev = [sPrev, s]; %#ok<AGROW>
              yPrev = [yPrev, y]; %#ok<AGROW>
              de    = (y'*y)/(y'*s);
            end
          end
          Hf = LbfgsProd(sPrev, yPrev, de);
        else
          sPrev = zeros(length(x), 0);
          yPrev = zeros(length(x), 0);
        end
    end
    
    % ------------ Solve subproblem for a search direction ------------
    
    if iter > 1 || evalHess
      quadF = @(z) smooth_quad(Hf, Df, x ,z);
      
      switch subMethod
        
        % Yuekai's implementation of SpaRSA
        case 'Sparsa'
          SparsaOptions = PNoptimset(SparsaOptions ,...
            'optTol', max( optTol, forcingTerm*opt ) ...
            );  
          
          [xProx, ~, spgOutput] = ...
            spg(quadF, nonsmoothF, x, SparsaOptions); 

          % ------------ Collect data from subproblem solve ------------
          
          subIters     = spgOutput.iters;
          subProxEvals = spgOutput.proxEvals;
          
          if debug
            subFlag = spgOutput.flag;
            subopt  = spgOutput.opt;
          end
        
        % TFOCS 
        case 'Tfocs'
          TfocsOpts.stopFcn = @(f,x) tfocs_stop( x, nonsmoothF,...
            max( optTol, forcingTerm*opt ) );

          [xProx, TfocsOut] = ...
            tfocs(quadF, [], nonsmoothF, x, TfocsOpts);
          
          subIters       = TfocsOut.niter;
          if isfield(TfocsOpts, 'countOps') && TfocsOpts.countOps
            subProxEvals = TfocsOut.counts(end,5);
          else
            subProxEvals = TfocsOut.niter;
          end
          
          if debug
            switch TfocsOut.status
              case 'Reached user''s supplied stopping criteria no. 1'
                subFlag = FLAG_OPT;
              case {'Step size tolerance reached (||dx||=0)',...
                    'Step size tolerance reached'           ,...
                    'Unexpectedly small stepsize'}
                subFlag = FLAG_TOLX;
              case 'Iteration limit reached'
                subFlag = FLAG_MAXITER;
              case 'Function/operator count limit reached'
                subFlag = FLAG_MAXFUNEVALS;
              otherwise
                subFlag = FLAG_OTHER;
            end
            subopt = TfocsOut.err(end);
          end
      end
      
      searchDir = xProx - x;
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
          descParam, xTol, maxfunEvals - funEvals); 
    else
      if iter > 1
        [x, f, Df, step, lineSearchFlag ,lineSearchIters] = ...
          LineSearch(x, searchDir, 1, f, h, Df'*searchDir, smoothF, nonsmoothF,...
            descParam, xTol, maxfunEvals - funEvals);
      else
        [x, f, Df, step, lineSearchFlag ,lineSearchIters] = ...
          CurvySearch(x, searchDir, max(min(1,1/norm(Df)),xTol), f, Df'*searchDir,...
            smoothF, nonsmoothF, descParam, xTol, maxfunEvals - funEvals); 
      end
    end
    
    % ------------ Select forcing term ------------
    
    if iter > 1 || evalHess
      [~, quadDf] = quadF(x);
      forcingTerm = min( 0.5, norm(Df-quadDf)/norm(Df) );
    end
    
    % ------------ Collect data for display and output ------------
    
    funEvals  =  funEvals + lineSearchIters ;
    proxEvals = proxEvals + lineSearchIters + subProxEvals;
    [~, xProx ] = nonsmoothF( x - Df ,1);
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
    elseif norm(x-xPrev,'inf')/max(1,norm(xPrev,'inf')) <= xTol 
      flag    = FLAG_TOLX;
      message = MSG_TOLX;
      loop    = 0;
    elseif abs(fPrev-f)/max(1,abs(fPrev)) <= funTol
      flag    = FLAG_TOLFUN;
      message = MSG_TOLFUN;
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
  
  clear global quadDf quadopt
  