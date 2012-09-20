function varargout = pnopt_backtrack( x, d, t, f_x, h_x, Df_Dx, smoothF, ...
  nonsmoothF, descParam, xtol, maxIter )
% pnopt_backtrack : Backtracking line search for step that satisfies a sufficient 
%   descent condition.
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $
% 
% --------------------Initialize--------------------

  % Set line search parameters
  be = 0.5;

  % Set termination flags
  FLAG_SUFFDESC = 1;
  FLAG_TOLX     = 2;
  FLAG_MAXFUNEV = 3;

  iter = 0;
  
  h_x1 = nonsmoothF( x + d );

  % --------------------Main Loop--------------------
  while 1
    iter = iter + 1;
    
    % Evaluate trial point and function value.
    x_t = x + d * t;
    if nargout > 6
      [ f_xt, Df_xt, Hf_xt ] = smoothF( x_t );
    else
      [ f_xt, Df_xt ] = smoothF( x_t );
    end
    h_xt = nonsmoothF( x_t );
    f_xt = f_xt + h_xt;
    
    % Check termination criteria
    desc = Df_Dx + h_x1 - h_x;
    if f_xt < f_x + descParam * t * desc         % Sufficient descent condition satisfied
      flag = FLAG_SUFFDESC;  
      break
    elseif t <= xtol            % Step length too small
      flag = FLAG_TOLX;
      break
    elseif iter >= maxIter      % Too many line search iterations
      flag = FLAG_MAXFUNEV;
      break
    end

    % Backtrack if objective value not well-defined
    if isnan( f_xt ) || isinf( f_xt ) || abs( f_xt - f_x-t * Df_Dx ) <= 1e-9
      t = be * t;
    % Safeguard quadratic interpolation
    else
      t_interp = - ( Df_Dx * t^2) / ( 2 * ( f_xt - f_x - t * Df_Dx ) );
      if 0.01 <= t_interp && t_interp <= 0.99*t 
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
  
  