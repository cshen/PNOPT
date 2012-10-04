function varargout = pnopt_curvtrack( x, d, t, f_old, df_x, smoothF, nonsmoothF, ...
  desc_param, xtol, maxIter )
% pnopt_curvtrack : Curve search for step that satisfies the Armijo condition
% 
%   $Revision: 0.1.1 $  $Date: 2012/06/30 $
% 
% ------------ Initialize ------------
  % Set line search parameters
  beta = 0.5;

  % Set termination flags
  FLAG_SUFFDESC = 1;
  FLAG_TOLX     = 2;
  FLAG_MAXFUNEV = 3;

  iter = 0;
  
  % ------------ Main Loop ------------
  while 1
    iter = iter + 1;
    
    % Evaluate trial point and function value.
    [ h_y, y ]  = nonsmoothF( x + t * d, t );
    if nargout > 6
      [ f_y, Df_y, D2f_y ] = smoothF( y );
    else
      [ f_y, Df_y ] = smoothF( y );
    end
    f_y = f_y + h_y;
    
    % Check termination criteria
    desc = 0.5 * norm( y - x ) ^2;
    if f_y < max( f_old ) + desc_param * t * desc    % Sufficient descent condition satisfied
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
    if isnan( f_y ) || isinf( f_y ) || abs( f_y - f_old(end) - t * df_x ) <= 1e-9
      t = beta * t;
    % Safeguard quadratic interpolation
    else
      t_interp = - ( df_x * t ^2) / ( 2 * ( f_y - f_old(end) - t * df_x ) );
      if 0.1 <= t_interp || t_interp <= 0.9*t 
        t = t_interp;
      else
        t = beta * t;
      end
    end
  end 
  
  if nargout > 6
    varargout = { y, f_y, Df_y, D2f_y, t, flag ,iter };
  else
    varargout = { y, f_y, Df_y, t, flag ,iter };
  end
  