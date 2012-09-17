function varargout = ...
  LineSearch(x, d, t, f, nonsmoothf, Dftd, smoothF, nonsmoothF, descParam, xtol, maxIter)
% LineSearch : Backtracking line search for step that satisfies a sufficient descent
%   condition.
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $
% 
% --------------------Initialize--------------------

  % Set line search parameters
  be = 0.5;

  % Set termination flags
  FLAG_SUFFDESCENT = 1;
  FLAG_TOLX        = 2;
  FLAG_MAXFUNEVALS = 3;

  iter = 0;
  
  Nonsmoothf1 = nonsmoothF(x+d);

  % --------------------Main Loop--------------------
  while 1
    iter = iter + 1;
    
    % Evaluate trial point and function value.
     xt = x+d*t;
    if nargout > 6
      [ft, Dft, Hft] = smoothF(xt);
    else
      [ft, Dft] = smoothF(xt);
    end
     ht = nonsmoothF(xt);
     ft = ft + ht;
    
    % Check termination criteria
    De = Dftd + Nonsmoothf1 - nonsmoothf;
    if ft < f + descParam*t*De         % Sufficient descent condition satisfied
      flag = FLAG_SUFFDESCENT;  
      break
    elseif t <= xtol            % Step length too small
      flag = FLAG_TOLX;
      break
    elseif iter >= maxIter      % Too many line search iterations
      flag = FLAG_MAXFUNEVALS;
      break
    end

    % Backtrack if objective value not well-defined
    if isnan(ft) || isinf(ft) || abs(ft-f-t*Dftd) <= 1e-9
      t = be*t;
    % Safeguard quadratic interpolation
    else
      tq = (-Dftd*t^2) / (2*(ft-f-t*Dftd));
      if 0.01 <= tq && tq <= 0.99*t 
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
  
  