function [ x, f_x, output ] = pnopt_PQN( smoothF, nonsmoothF, x, options )
% pnopt_PQN : Proximal quasi-Newton methods
% 
%   $Revision: 0.8.0 $  $Date: 2012/12/01 $
% 
  REVISION = '$Revision: 0.8.0$';
  DATE     = '$Date: Dec. 01, 2012$';
  REVISION = REVISION(11:end-1);
  DATE     = DATE(8:end-1);
  
% ============ Process options ============
  
  sparsa_options = pnopt_optimset(...
    'display'  , 0    ,...
    'maxfunEv' , 5000 ,...
    'maxIter'  , 500   ...
    );
  
  tfocsOpts = struct(...
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
  method         = options.method;
  subProb_solver = options.subProb_solver;
  ftol           = options.ftol;
  optim_tol      = options.optim_tol;
  xtol           = options.xtol;
  
  switch method
    case 'bfgs'
      
    case 'Lbfgs'
      Lbfgs_mem = options.Lbfgs_mem;
  end
  
% ------------ Set subproblem solver options ------------
  
  switch subProb_solver
    case 'sparsa'
      if isfield( options, 'sparsa_options' ) && ~isempty( options.sparsa_options )
        sparsa_options = pnopt_optimset( sparsa_options, options.sparsa_options );
      end
    case 'tfocs'
      if isfield( options, 'tfocsOpts' ) && ~isempty( options.tfocsOpts )
        tfocsOpts = merge_struct( tfocsOpts, options.tfocsOpts );
      end
           
      if debug
        tfocsOpts.countOps = 1;
        tfocsOpts.errFcn   = @ ( f, x ) tfocs_err();
      end
  end
  
% ============ Initialize variables ============
  
  pnopt_flags
  
  iter         = 0; 
  loop         = 1;
  forcing_term = 0.1;
  
  Trace.f_x    = zeros( maxIter + 1, 1 );
  Trace.funEv  = zeros( maxIter + 1, 1 );
  Trace.proxEv = zeros( maxIter + 1, 1 );
  Trace.optim  = zeros( maxIter + 1, 1 );
  
  if debug
    Trace.forcing_term    = zeros( maxIter, 1 );
    Trace.subProb_iters   = zeros( maxIter, 1 );
    Trace.subProb_optim   = zeros( maxIter, 1 );
  end
  
  if display > 0  
    fprintf( ' %s\n', repmat( '=', 1, 64 ) );
    fprintf( '                  PNOPT v.%s (%s)\n', REVISION, DATE );
    fprintf( ' %s\n', repmat( '=', 1, 64 ) );
    fprintf( ' %4s   %6s  %6s  %12s  %12s  %12s \n',...
      '','Fun.', 'Prox', 'Step len.', 'Obj. val.', 'Optim.' );
    fprintf( ' %s\n', repmat( '-', 1, 64 ) );
  end
  
% ------------ Evaluate objective function at starting x ------------
  
  [ g_x, Dg_x ] = smoothF( x );
    h_x         = nonsmoothF( x );
    f_x         = g_x + h_x;
  
% ------------ Start collecting data for display and output ------------
  
    funEv       = 1;
    proxEv      = 0;
  [ ~, x_prox ] = nonsmoothF( x - Dg_x, 1 );
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
    
  % ------------ Update Hessian approximation ------------
    
    switch method
      
      % BFGS method
      case 'bfgs'
        if iter > 1
          s    =  x - x_old;
          y    = Dg_x - Dg_old;
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
          y = Dg_x - Dg_old;
          if y'*s > 1e-9
            if size( sPrev, 2 ) > Lbfgs_mem
              sPrev = [ sPrev(:,2:Lbfgs_mem), s ];
              yPrev = [ yPrev(:,2:Lbfgs_mem), y ];
              hDiag = ( y' * y ) / ( y' * s );
            else
              sPrev = [ sPrev, s ]; %#ok<AGROW>
              yPrev = [ yPrev, y ]; %#ok<AGROW>
              hDiag = ( y' * y ) / ( y' * s );
            end
          end
          Hf_x = pnopt_Lbfgs_prod( sPrev, yPrev, hDiag );
        else 
          sPrev = zeros( length(x), 0 );
          yPrev = zeros( length(x), 0 );
          hDiag = 1;
        end
    end
    
  % ------------ Solve subproblem for a search direction ------------
    
    if iter > 1 
      quadF = @(z) pnopt_quad( Hf_x, Dg_x, f_x, z - x );
      
      switch subProb_solver
        
        % SpaRSA
        case 'sparsa'
          sparsa_options = pnopt_optimset( sparsa_options ,...
            'optim_tol', max( 0.1 * optim_tol, forcing_term * optim ) ...
            );  
          
          [ x_prox, ~, sparsa_out ] = ...
            pnopt_sparsa( quadF, nonsmoothF, x, sparsa_options ); 

        % ------------ Collect data from subproblem solve ------------
          
          subProb_iters  = sparsa_out.iters;
          subProb_proxEv = sparsa_out.proxEv;
          
          if debug
            subProb_flag  = sparsa_out.flag; %#ok<NASGU>
            subProb_optim = sparsa_out.optim;
          end
        
        % TFOCS 
        case 'tfocs'
          tfocsOpts.stopFcn = @(f, x) tfocs_stop( x, nonsmoothF,...
            max( 0.1 * optim_tol, forcing_term * optim ) );

          [ x_prox, tfocsOut ] = ...
            tfocs( quadF, [], nonsmoothF, x, tfocsOpts );
        
          subProb_iters = tfocsOut.niter;
          if isfield( tfocsOpts, 'countOps' ) && tfocsOpts.countOps
            subProb_proxEv = tfocsOut.counts(end,5);
          else
            subProb_proxEv = tfocsOut.niter;
          end
          
          if debug
            subProb_optim = tfocsOut.err(end);
          end
      end
      
      Dx = x_prox - x;
    else
      subProb_iters  = 0;
      subProb_proxEv = 0;
      
      if debug 
        subProb_optim = 0;
      end
      
      Dx = - Dg_x;
    end
    
  % ------------ Conduct line search ------------
    
    x_old  = x;
    f_old  = f_x;
    Dg_old = Dg_x;
    
    if iter > 1
      [ x, f_x, Dg_x, step, backtrack_flag, backtrack_iters ] = ...
        pnopt_backtrack( x, Dx, 1, f_x, h_x, Dg_x' * Dx, smoothF, nonsmoothF, ...
          desc_param, xtol, maxfunEv - funEv ); %#ok<ASGLU>
    else
      [ x, f_x, Dg_x, step, backtrack_flag, backtrack_iters ] = ...
        pnopt_curvtrack( x, Dx, max( min( 1, 1 / norm( Dg_x ) ), xtol ), f_x, ...
          Dg_x'*Dx, smoothF, nonsmoothF, desc_param, xtol, maxfunEv - funEv );  %#ok<ASGLU>
    end
    
  % ------------ Select safeguarded forcing term ------------
    
    if iter > 1 
      [ q_x, subProb_Dg_x ] = quadF(x);  %#ok<ASGLU>
        forcing_term        = min( 0.1 , norm( Dg_x - subProb_Dg_x ) / norm( Dg_old ) );
    end
    
  % ------------ Collect data for display and output ------------
    
      funEv       =  funEv + backtrack_iters ;
      proxEv      = proxEv + backtrack_iters + subProb_proxEv;
    [ ~, x_prox ] = nonsmoothF( x - Dg_x, 1 );
      optim       = norm( x_prox - x, 'inf' );
    
    Trace.f_x(iter+1)    = f_x;
    Trace.funEv(iter+1)  = funEv;
    Trace.proxEv(iter+1) = proxEv;
    Trace.optim(iter+1)  = optim;
    
    if debug
      Trace.forcing_term(iter)    = forcing_term;
      Trace.backtrack_iters(iter) = backtrack_iters;
      Trace.subProb_iters(iter)   = subProb_iters;
      Trace.subProb_optim(iter)   = subProb_optim;
    end
    
    if display > 0 && mod( iter, display ) == 0
      fprintf( ' %4d | %6d  %6d  %12.4e  %12.4e  %12.4e\n', ...  
        iter, funEv, proxEv, step, f_x, optim );
    end
    
    pnopt_stop
    
  end
  
% ============ Clean up and exit ============
  
  Trace.f_x    = Trace.f_x(1:iter+1);
  Trace.funEv  = Trace.funEv(1:iter+1);
  Trace.proxEv = Trace.proxEv(1:iter+1);
  Trace.optim  = Trace.optim(1:iter+1);
  
  if debug
    Trace.forcing_term    = Trace.forcing_term(1:iter);
    Trace.backtrack_iters = Trace.backtrack_iters(1:iter);
    Trace.subProb_iters   = Trace.subProb_iters(1:iter);
    Trace.subProb_optim   = Trace.subProb_optim(1:iter);
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
  
  clear global subProb_Dg_y subProb_optim
    