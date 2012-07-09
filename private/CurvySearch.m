function varargout = ...
  CurvySearch(x, d, t, f, gtd, smoothF, nonsmoothF, descCond, TolX, maxIter)
% CurvySearch : Curve search for step that satisfies the Armijo condition
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $
% 
% ------------ Initialize ------------
  % Set line search parameters
  be = 0.5;

  % Set termination flags
  FLAG_SUFFDESCENT = 1;
  FLAG_TOLX        = 2;
  FLAG_MAXFUNEVALS = 3;

  iter = 0;
  
  % ------------ Main Loop ------------
  while 1
    iter = iter + 1;
    
    % Evaluate trial point and function value.
    [ht, xt]  = nonsmoothF(x+t*d, t);
    if nargout > 6
      [ft, Dft, Hft] = smoothF(xt);
    else
      [ft, Dft] = smoothF(xt);
    end
    ft = ft + ht;
    
    % Check termination criteria
    De = 0.5*norm(xt-x)^2;
    if ft < max(f) + descCond*t*De    % Sufficient descent condition satisfied
      flag = FLAG_SUFFDESCENT;  
      break
    elseif t <= TolX            % Step length too small
      flag = FLAG_TOLX;
      break
    elseif iter >= maxIter      % Too many line search iterations
      flag = FLAG_MAXFUNEVALS;
      break
    end

    % Backtrack if objective value not well-defined of function seems linear
    if isnan(ft) || isinf(ft) || abs(ft-f(end)-t*gtd) <= 1e-9
      t = be*t;
    % Safeguard quadratic interpolation
    else
      tq = (-gtd*t^2) / (2*(ft-f(end)-t*gtd));
      if 0.1 <= tq || tq <= 0.9*t 
        t = tq;
      else
        t = be*t;
      end
    end
  end 
  
  if nargout > 6
    varargout = {xt, ft, Dft, Hft, t, flag ,iter};
  else
    varargout = {xt, ft, Dft, t, flag ,iter};
  end
  