function varargout = pnopt_curvtrack( x, d, t, f_old, Df_Dx, smoothF, nonsmoothF, ...
  descParam, xtol, maxIter )
% pnopt_curvtrack : Curve search for step that satisfies the Armijo condition
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $
% 
% ------------ Initialize ------------
  % Set line search parameters
  be = 0.5;

  % Set termination flags
  FLAG_SUFFDESC = 1;
  FLAG_TOLX     = 2;
  FLAG_MAXFUNEV = 3;

  iter = 0;
  
  % ------------ Main Loop ------------
  while 1
    iter = iter + 1;
    
    % Evaluate trial point and function value.
    [ h_xt, x_t ]  = nonsmoothF( x + t * d, t );
    if nargout > 6
      [ f_xt, Df_xt, Hf_xt ] = smoothF( x_t );
    else
      [ f_xt, Df_xt ] = smoothF( x_t );
    end
    f_xt = f_xt + h_xt;
    
    % Check termination criteria
    desc = 0.5 * norm( x_t - x ) ^2;
    if f_xt < max( f_old ) + descParam * t * desc    % Sufficient descent condition satisfied
      flag = FLAG_SUFFDESC;  
      break
    elseif t <= xtol            % Step length too small
      flag = FLAG_TOLX;
      break
    elseif iter >= maxIter      % Too many line search iterations
      flag = FLAG_MAXFUNEV;
      break
    end

    % Backtrack if objective value not well-defined of function seems linear
    if isnan( f_xt ) || isinf( f_xt ) || abs( f_xt - f_old(end) - t * Df_Dx ) <= 1e-9
      t = be * t;
    % Safeguard quadratic interpolation
    else
      t_interp = - ( Df_Dx * t ^2) / ( 2 * ( f_xt - f_old(end) - t * Df_Dx ) );
      if 0.1 <= t_interp || t_interp <= 0.9*t 
        t = t_interp;
      else
        t = be * t;
      end
    end
  end 
  
  if nargout > 6
    varargout = { x_t, f_xt, Df_xt, Hf_xt, t, flag ,iter };
  else
    varargout = { x_t, f_xt, Df_xt, t, flag ,iter };
  end
  