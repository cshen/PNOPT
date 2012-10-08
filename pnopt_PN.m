function [ x, f_x, output ] = pnopt_PN( smoothF, nonsmoothF, x, options )
% pnopt_PN : Proximal Newton method
% 
%   $Revision: 0.6.4 $  $Date: 2012/09/30 $
% 
  REVISION = '$Revision: 0.5.1$';
  DATE     = '$Date: Sep. 15, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Process options ============
  
  sparsa_options = pnopt_optimset(...
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
  
  debug          = options.debug;
  desc_param     = options.desc_param;
  display        = options.display;
  maxfunEv       = options.maxfunEv;
  maxIter        = options.maxIter;
  subprob_solver = options.subprob_solver;
  ftol           = options.ftol;
  optim_tol      = options.optim_tol;
  xtol           = options.xtol;
  
% ------------ Set subproblem solver options ------------
  
  switch subprob_solver
    case 'sparsa'
      if isfield( options, 'sparsaOpts' ) && ~isempty( options.sparsaOpts )
        sparsa_options = pnopt_optimset( sparsa_options, options.sparsaOpts );
      end
    case 'tfocs'
      if isfield( options, 'TfocsOpts' ) && ~isempty( options.TfocsOpts )
        TfocsOpts = merge_struct( TfocsOpts, options.TfocsOpts );
      end
           
      if debug
        TfocsOpts.countOps = 1;
        TfocsOpts.errFcn   = @(f, x) tfocs_err();
      end
  end
  
% ============ Initialize variables ============
  
  FLAG_OPTIM   = 1;
  FLAG_XTOL    = 2;
  FLAG_FTOL    = 3;
  FLAG_MAXITER = 4;
  FLAG_MAXFEV  = 5;
  FLAG_OTHER   = 6;
  
  MESSAGE_OPTIM   = 'Optimality below optim_tol.';
  MESSAGE_XTOL    = 'Relative change in x below xtol.';
  MESSAGE_FTOL    = 'Relative change in function value below funTol.';
  MESSAGE_MAXITER = 'Max number of iterations reached.';
  MESSAGE_MAXFEV  = 'Max number of function evaluations reached.';
  
  iter         = 0; 
  loop         = 1;
  forcing_term = 0.5;
  
  Trace.f_x    = zeros( maxIter+1, 1 );
  Trace.funEv  = zeros( maxIter+1, 1 );
  Trace.proxEv = zeros( maxIter+1, 1 );
  Trace.optim  = zeros( maxIter+1, 1 );
  
  if debug
    Trace.forcing_term    = zeros( maxIter, 1 );
    Trace.normDx          = zeros( maxIter, 1 );
    Trace.backtrack_flag  = zeros( maxIter, 1 );
    Trace.backtrack_iters = zeros( maxIter, 1 );
    Trace.subprob_flags   = zeros( maxIter, 1 );
    Trace.subprob_iters   = zeros( maxIter, 1 );
    Trace.subprob_optim   = zeros( maxIter, 1 );
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
  
  [ f_x, Df_x, D2f_x ] = smoothF( x );
    h_x                = nonsmoothF( x );
    f_x                = f_x + h_x;
  
% ------------ Start collecting data for display and output ------------
  
    funEv       = 1;
    proxEv      = 0;
  [ ~, x_prox ] = nonsmoothF( x - Df_x, 1 );
    optim       = norm( x_prox - x, 'inf' );
  
  Trace.f_x(1)    = f_x;
  Trace.funEv(1)  = funEv;
  Trace.proxEv(1) = proxEv;
  Trace.optim(1)  = optim; 
  
  if display > 0
    fprintf( ' %4d | %6d  %6d  %12s  %12.4e  %12.4e\n', ...
      iter, funEv, proxEv, '', f_x, optim );
  end
  
% ------------ Check if starting x is optimal ------------
  
  if optim <= optim_tol
    flag    = FLAG_OPTIM;
    message = MESSAGE_OPTIM;
    loop    = 0;
  end

% ============ Main Loop ============
  
  while loop
    iter = iter + 1; 
    
  % ------------ Solve subproblem for a search direction ------------
    
    quadF = @(z) pnopt_quad( D2f_x, Df_x, f_x, z - x );

    switch subprob_solver

      % SpaRSA
      case 'sparsa'
        sparsa_options = pnopt_optimset( sparsa_options ,...
          'optim_tol', max( 0.5 * optim_tol, forcing_term * optim ) ...
          );  

        [ x_prox, ~, sparsa_out ] = ...
          pnopt_sparsa( quadF, nonsmoothF, x, sparsa_options ); 

      % ------------ Collect data from subproblem solve ------------

        subprob_iters  = sparsa_out.iters;
        subprob_proxEv = sparsa_out.proxEv;

        if debug
          subprob_flag  = sparsa_out.flag;
          subprob_optim = sparsa_out.optim;
        end

      % TFOCS 
      case 'tfocs'
        TfocsOpts.stopFcn = @(f, x) tfocs_stop( x, nonsmoothF,...
          max( 0.5 * optim_tol, forcing_term * optim ) );

        [ x_prox, tfocsOut ] = ...
          tfocs( quadF, [], nonsmoothF, x, TfocsOpts );

        subprob_iters = tfocsOut.niter;
        if isfield( TfocsOpts, 'countOps' ) && TfocsOpts.countOps
          subprob_proxEv = tfocsOut.counts(end,5);
        else
          subprob_proxEv = tfocsOut.niter;
        end

        if debug
          switch tfocsOut.status
            case 'Reached user''s supplied stopping criteria no. 1'
              subprob_flag = FLAG_OPTIM;
            case { 'Step size tolerance reached (||dx||=0)', ...
                   'Step size tolerance reached'           , ...
                   'Unexpectedly small stepsize' }
              subprob_flag = FLAG_XTOL;
            case 'Iteration limit reached'
              subprob_flag = FLAG_MAXITER;
            case 'Function/operator count limit reached'
              subprob_flag = FLAG_MAXFEV;
            otherwise
              subprob_flag = FLAG_OTHER;
          end
          subprob_optim = tfocsOut.err(end);
        end
    end
      
    Dx = x_prox - x;
    
  % ------------ Conduct line search ------------
    
    x_old  = x;
    f_old  = f_x;
    
    [ x, f_x, Df_x, D2f_x, step, backtrack_flag, backtrack_iters ] = ...
      pnopt_backtrack( x, Dx, 1, f_x, h_x, Df_x' * Dx, smoothF, nonsmoothF, ...
        desc_param, xtol, maxfunEv - funEv ); 
    
  % ------------ Select safeguarded forcing term ------------
    
    [ quadf, subprob_Df_x ] = quadF(x);  %#ok<ASGLU>
      forcing_term       = min( 0.5, norm( Df_x - subprob_Df_x ) / norm( Df_x ) );
    
  % ------------ Collect data for display and output ------------
    
      funEv       =  funEv + backtrack_iters ;
      proxEv      = proxEv + backtrack_iters + subprob_proxEv;
    [ ~, x_prox ] = nonsmoothF( x - Df_x, 1 );
      optim         = norm( x_prox - x, 'inf' );
    
    Trace.f_x(iter+1)    = f_x;
    Trace.funEv(iter+1)  = funEv;
    Trace.proxEv(iter+1) = proxEv;
    Trace.optim(iter+1)  = optim;
    
    if debug
      Trace.forc_term(iter)       = forcing_term;
      Trace.normDx(iter)          = norm(Dx);
      Trace.backtrack_flag(iter)  = backtrack_flag;
      Trace.backtrack_iters(iter) = backtrack_iters;
      Trace.subprob_flags(iter)   = subprob_flag;
      Trace.subprob_iters(iter)   = subprob_iters;
      Trace.subprob_optim(iter)   = subprob_optim;
    end
    
    if display > 0 && mod( iter, display ) == 0
      fprintf( ' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n', ...  
        iter, funEv, proxEv, step, f_x, optim );
    end
    
  % ------------ Check stopping criteria ------------
    
    if optim <= optim_tol
      flag    = FLAG_OPTIM;
      message = MESSAGE_OPTIM;
      loop    = 0;
    elseif norm( x - x_old, 'inf' ) / max( 1, norm( x_old, 'inf' ) ) <= xtol 
      flag    = FLAG_XTOL;
      message = MESSAGE_XTOL;
      loop    = 0;
    elseif abs( f_old - f_x ) / max( 1, abs( f_old ) ) <= ftol
      flag    = FLAG_FTOL  ;
      message = MESSAGE_FTOL  ;
      loop    = 0;
    elseif iter >= maxIter 
      flag    = FLAG_MAXITER;
      message = MESSAGE_MAXITER;
      loop    = 0;
    elseif funEv >= maxfunEv
      flag    = FLAG_MAXFEV;
      message = MESSAGE_MAXFEV;
      loop    = 0;
    end
  end
  
% ============ Cleanup and exit ============
  
  Trace.f_x    = Trace.f_x(1:iter+1);
  Trace.funEv  = Trace.funEv(1:iter+1);
  Trace.proxEv = Trace.proxEv(1:iter+1);
  Trace.optim  = Trace.optim(1:iter+1);
  
  if debug
    Trace.forcing_term    = Trace.forc_term(1:iter);
    Trace.normDx          = Trace.normDx(1:iter);
    Trace.backtrack_flag  = Trace.backtrack_flag(1:iter);
    Trace.backtrack_iters = Trace.backtrack_iters(1:iter);
    Trace.subprob_flags   = Trace.subprob_flags(1:iter);
    Trace.subprob_iters   = Trace.subprob_iters(1:iter);
    Trace.subprob_optim   = Trace.subprob_optim(1:iter);
  end
  
  if display > 0 && mod( iter, display ) > 0
    fprintf( ' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n', ...
      iter, funEv, proxEv, step, f_x, optim );
  end
  
  output = struct( ...    
    'flag'    , flag    ,...
    'funEv'   , funEv   ,...
    'iters'   , iter    ,...
    'optim'   , optim   ,...
    'options' , options ,...
    'proxEv'  , proxEv  ,...
    'Trace'   , Trace    ...
    );
  
  if display > 0
    fprintf( ' %s\n', repmat( '-', 1, 64 ) );
    fprintf( ' %s\n', message )
    fprintf( ' %s\n', repmat( '-', 1, 64 ) );
  end
  
  clear global Dq_y subprob_optim
  