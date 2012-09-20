function [x, f_x, output] = pnopt_PQN(smoothF, nonsmoothF, x, options)
% pnopt_PQN : Proximal quasi-Newton methods
% 
%   $Revision: 0.5.1 $  $Date: 2012/09/15 $
% 
  REVISION = '$Revision: 0.5.1$';
  DATE     = '$Date: Sep. 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Process options ============
  
  SparsaOpts = pnopt_optimset(...
    'display'  , 0    ,...
    'maxfunEv' , 5000 ,...
    'maxIter'  , 500   ...
    );
  
  TfocsOpts = struct(...
    'alg'        , 'N83' ,...
    'maxIts'     , 500   ,...
    'printEvery' , 0     ,...
    'restart'    , -Inf   ...
    );
  
  debug       = options.debug;
  desc_param  = options.desc_param;
  display     = options.display;
  maxfunEv    = options.maxfunEv;
  maxIter     = options.maxIter;
  method      = options.method;
  quad_solver = options.quad_solver;
  funTol      = options.funTol;
  optTol      = options.optTol;
  xtol        = options.xtol;
  
  switch method
    case 'bfgs'
      
    case 'Lbfgs'
      Lbfgs_mem = options.Lbfgs_mem;
  end
  
% ------------ Set subproblem quad_solver options ------------
  
  switch quad_solver
    case 'Sparsa'
      if isfield( options, 'SparsaOpts' ) && ~isempty( options.SparsaOpts )
        SparsaOpts = pnopt_optimset( SparsaOpts, options.SparsaOpts );
      end
    case 'Tfocs'
      if isfield( options, 'TfocsOpts' ) && ~isempty( options.TfocsOpts )
        TfocsOpts = merge_struct( TfocsOpts, options.TfocsOpts );
      end
           
      if debug
        TfocsOpts.countOps = 1;
        TfocsOpts.errFcn   = @(f, x) tfocs_err();
      end
  end
  
% ============ Initialize variables ============
  
  FLAG_OPT     = 1;
  FLAG_XTOL    = 2;
  FLAG_FUNTOL  = 3;
  FLAG_MAXITER = 4;
  FLAG_MAXFEV  = 5;
  FLAG_OTHER   = 6;
  
  MSG_OPT     = 'Optimality below optTol.';
  MSG_XTOL    = 'Relative change in x below xtol.';
  MSG_FUNTOL  = 'Relative change in function value below funTol.';
  MSG_MAXITER = 'Max number of iterations reached.';
  MSG_MAXFEV  = 'Max number of function evaluations reached.';
  
  iter         = 0; 
  loop         = 1;
  forcing_term = 0.5;
  
  Trace.f_x    = zeros( maxIter + 1, 1 );
  Trace.funEv  = zeros( maxIter + 1, 1 );
  Trace.proxEv = zeros( maxIter + 1, 1 );
  Trace.opt    = zeros( maxIter + 1, 1 );
  
  if debug
    Trace.forc_term       = zeros( maxIter, 1 );
    Trace.normDx          = zeros( maxIter, 1 );
    Trace.backtrack_flag  = zeros( maxIter, 1 );
    Trace.backtrack_iters = zeros( maxIter, 1 );
    Trace.quad_flags      = zeros( maxIter, 1 );
    Trace.quad_iters      = zeros( maxIter, 1 );
    Trace.quad_opt        = zeros( maxIter, 1 );
  end
  
  if display > 0  
    fprintf( ' %s\n', repmat( '=', 1, 64 ) );
    fprintf( '                  PNOPT v.%s (%s)\n', REVISION, DATE );
    fprintf( ' %s\n', repmat( '=', 1, 64 ) );
    fprintf( ' %4s   %6s  %6s  %12s  %12s  %12s \n',...
      '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'Optimality' );
    fprintf( ' %s\n', repmat( '-', 1, 64 ) );
  end
  
% ------------ Evaluate objective function at starting x ------------
  
  [ f_x, Df_x ] = smoothF( x );
    h_x         = nonsmoothF( x );
    f_x         = f_x + h_x;
  
% ------------ Start collecting data for display and output ------------
  
    funEv       = 1;
    proxEv      = 0;
  [ ~, x_prox ] = nonsmoothF( x - Df_x, 1 );
    opt         = norm( x_prox - x, 'inf' );
  
  Trace.f_x(1)    = f_x;
  Trace.funEv(1)  = funEv;
  Trace.proxEv(1) = proxEv;
  Trace.opt(1)    = opt; 
  
  if display > 0
    fprintf( ' %4d | %6d  %6d  %12s  %12.4e  %12.4e\n', ...
      iter, funEv, proxEv, '', f_x, opt );
  end
  
% ------------ Check if starting x is optimal ------------
  
  if opt <= optTol
    flag    = FLAG_OPT;
    message = MSG_OPT;
    loop    = 0;
  end

% ============ Main Loop ============
  
  while loop
    iter = iter + 1; 
    
  % ------------ Update Hessian approximation ------------
    
    switch method
      
      % BFGS method
      case 'bfgs'
        if iter > 1
          s    =  x - x_old;
          y    = Df_x - Df_old;
          qty1 = cholB' * ( cholB * s );
          if s'*y > 1e-9
            cholB = cholupdate( cholupdate( cholB, y / sqrt( y' *s ) ), qty1 / ...
              sqrt( s' * qty1 ), '-' );
          end
          Hf_x = @(x) cholB' * ( cholB * x );
        else
          cholB = eye( length( x ) );
        end

      % Limited-memory BFGS method
      case 'Lbfgs'
        if iter > 1
          s =  x - x_old;
          y = Df_x - Df_old;
          if y'*s > 1e-9
            if size( sPrev, 2 ) > Lbfgs_mem
              sPrev = [ sPrev(:,2:Lbfgs_mem), s ];
              yPrev = [ yPrev(:,2:Lbfgs_mem), y ];
              de    = ( y' * y ) / ( y' * s );
            else
              sPrev = [ sPrev, s ]; %#ok<AGROW>
              yPrev = [ yPrev, y ]; %#ok<AGROW>
              de    = ( y' * y ) / ( y' * s );
            end
          end
          Hf_x = pnopt_Lbfgs_prod( sPrev, yPrev, de );
        else 
          sPrev = zeros( length(x), 0 );
          yPrev = zeros( length(x), 0 );
          de    = 1;
        end
    end
    
  % ------------ Solve subproblem for a search direction ------------
    
    if iter > 1 
      quadF = @(z) smooth_quad( Hf_x, Df_x, f_x, z - x );
      
      switch quad_solver
        
        % SpaRSA
        case 'Sparsa'
          SparsaOpts = pnopt_optimset( SparsaOpts ,...
            'optTol', max( 0.5 * optTol, forcing_term * opt ) ...
            );  
          
          [ x_prox, ~, SparsaOut ] = ...
            pnopt_parsa( quadF, nonsmoothF, x, SparsaOpts ); 

        % ------------ Collect data from subproblem solve ------------
          
          quad_iters  = SparsaOut.iters;
          quad_proxEv = SparsaOut.proxEv;
          
          if debug
            quad_flag = SparsaOut.flag;
            quad_opt  = SparsaOut.opt;
          end
        
        % TFOCS 
        case 'Tfocs'
          TfocsOpts.stopFcn = @(f, x) tfocs_stop( x, nonsmoothF,...
            max( 0.5*optTol, forcing_term*opt ) );

          [ x_prox, TfocsOut ] = ...
            tfocs( quadF, [], nonsmoothF, x, TfocsOpts );
        
          quad_iters = TfocsOut.niter;
          if isfield( TfocsOpts, 'countOps' ) && TfocsOpts.countOps
            quad_proxEv = TfocsOut.counts(end,5);
          else
            quad_proxEv = TfocsOut.niter;
          end
          
          if debug
            switch TfocsOut.status
              case 'Reached user''s supplied stopping criteria no. 1'
                quad_flag = FLAG_OPT;
              case { 'Step size tolerance reached (||dx||=0)', ...
                     'Step size tolerance reached'           , ...
                     'Unexpectedly small stepsize' }
                quad_flag = FLAG_XTOL;
              case 'Iteration limit reached'
                quad_flag = FLAG_MAXITER;
              case 'Function/operator count limit reached'
                quad_flag = FLAG_MAXFEV;
              otherwise
                quad_flag = FLAG_OTHER;
            end
            quad_opt = TfocsOut.err(end);
          end
      end
      
      Dx = x_prox - x;
    else
      quad_iters  = 0;
      quad_proxEv = 0;
      
      if debug 
        quad_flag = 0;
        quad_opt  = 0;
      end
      
      Dx = - Df_x;
    end
    
  % ------------ Conduct line search ------------
    
    x_old  = x;
    f_old  = f_x;
    Df_old = Df_x;
    
    if iter > 1
      [ x, f_x, Df_x, step, backtrack_flag,backtrack_iters ] = ...
        pnopt_backtrack( x, Dx, 1, f_x, h_x, Df_x' * Dx, smoothF, nonsmoothF, ...
          desc_param, xtol, maxfunEv - funEv );
    else
      [ x, f_x, Df_x, step, backtrack_flag,backtrack_iters ] = ...
        pnopt_curvtrack( x, Dx, max( min( 1, 1 / norm( Df_x ) ), xtol ), f_x, ...
          Df_x'*Dx, smoothF, nonsmoothF, desc_param, xtol, maxfunEv - funEv ); 
    end
    
  % ------------ Select safeguarded forcing term ------------
    
    if iter > 1 
      [ q_x, quad_Df_x ] = quadF(x);  %#ok<ASGLU>
        forcing_term     = min( 0.5, norm( Df_x - quad_Df_x ) / norm( Df_x ) );
    end
    
  % ------------ Collect data for display and output ------------
    
      funEv       =  funEv + backtrack_iters ;
      proxEv      = proxEv + backtrack_iters + quad_proxEv;
    [ ~, x_prox ] = nonsmoothF( x - Df_x, 1 );
      opt         = norm( x_prox - x, 'inf' );
    
    Trace.f_x(iter+1)    = f_x;
    Trace.funEv(iter+1)  = funEv;
    Trace.proxEv(iter+1) = proxEv;
    Trace.opt(iter+1)    = opt;
    
    if debug
      Trace.forc_term(iter)       = forcing_term;
      Trace.normDx(iter)          = norm(Dx);
      Trace.backtrack_flag(iter)  = backtrack_flag;
      Trace.backtrack_iters(iter) = backtrack_iters;
      Trace.quad_flags(iter)      = quad_flag;
      Trace.quad_iters(iter)      = quad_iters;
      Trace.quad_opt(iter)        = quad_opt;
    end
    
    if display > 0 && mod( iter, display ) == 0
      fprintf( ' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n', ...  
        iter, funEv, proxEv, step, f_x, opt );
    end
    
  % ------------ Check stopping criteria ------------
    
    if opt <= optTol
      flag    = FLAG_OPT;
      message = MSG_OPT;
      loop    = 0;
    elseif norm( x - x_old, 'inf' ) / max( 1, norm( x_old, 'inf' ) ) <= xtol 
      flag    = FLAG_XTOL;
      message = MSG_XTOL;
      loop    = 0;
    elseif abs( f_old - f_x ) / max( 1, abs( f_old ) ) <= funTol
      flag    = FLAG_FUNTOL;
      message = MSG_FUNTOL;
      loop    = 0;
    elseif iter >= maxIter 
      flag    = FLAG_MAXITER;
      message = MSG_MAXITER;
      loop    = 0;
    elseif funEv >= maxfunEv
      flag    = FLAG_MAXFEV;
      message = MSG_MAXFEV;
      loop    = 0;
    end
  end
  
% ============ Cleanup and exit ============
  
  Trace.f_x    = Trace.f_x(1:iter+1);
  Trace.funEv  = Trace.funEv(1:iter+1);
  Trace.proxEv = Trace.proxEv(1:iter+1);
  Trace.opt    = Trace.opt(1:iter+1);
  
  if debug
    Trace.forc_term       = Trace.forc_term(1:iter);
    Trace.normDx          = Trace.normDx(1:iter);
    Trace.backtrack_flag  = Trace.backtrack_flag(1:iter);
    Trace.backtrack_iters = Trace.backtrack_iters(1:iter);
    Trace.quad_flags      = Trace.quad_flags(1:iter);
    Trace.quad_iters      = Trace.quad_iters(1:iter);
    Trace.quad_opt        = Trace.quad_opt(1:iter);
  end
  
  if display > 0 && mod( iter, display ) > 0
    fprintf( ' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n', ...
      iter, funEv, proxEv, step, f_x, opt );
  end
  
  output = struct( ...
    'flag'    , flag      ,...
    'funEv'   , funEv  ,...
    'iters'   , iter      ,...
    'opt'     , opt       ,...
    'options' , options   ,...
    'proxEv'  , proxEv ,...
    'Trace'   , Trace      ...
    );
  
  if display > 0
    fprintf( ' %s\n', repmat( '-', 1, 64 ) );
    fprintf( ' %s\n', message )
    fprintf( ' %s\n', repmat( '-', 1, 64 ) );
  end
  
  clear global quad_Df_x quad_opt
  