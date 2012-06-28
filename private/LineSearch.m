function [xt, ft, ht, Dft, t, Flag ,iter] = LineSearch ...
  (x, d, t, f, nonsmoothf, gtd, smoothF, nonsmoothF, TolX, maxIter)
% LineSearch : Line search for step that satisfies a sufficient descent
%   condition.
% 
%   $Revision: 0.1.0 $  $Date: 2012/05/30 $

% --------------------Initialize--------------------
  % Set line search parameters
  al = 0.0001;
  be = 0.5;

  % Set termination flags
  FLAG_SUFFDESCENT = 1;
  FLAG_TOLX        = 2;
  FLAG_MAXFUNEVALS = 3;

  iter = 0;
  
  Nonsmoothf1 = nonsmoothF(x+d);

  % --------------------Main Loop--------------------
  while 1
    % Evaluate trial point and function value.
     xt = x+d*t;
    [ft, Dft] = smoothF(xt);
     ht       = nonsmoothF(xt);
     ft       = ft + ht;
     
    iter = iter + 1;
    
    % Check termination criteria
    De = gtd + Nonsmoothf1 - nonsmoothf;
    if ft < f + al*t*De             % Sufficient descent condition satisfied
      Flag = FLAG_SUFFDESCENT;  
      break
    elseif t <= TolX                % Step length too small
      Flag = FLAG_TOLX;
      break
    elseif iter >= maxIter  % Too many linesearch iterations.
      Flag = FLAG_MAXFUNEVALS;
      break
    end

    % Backtrack if objective value not well-defined
    if isnan(ft) || isinf(ft) || abs(ft-f-t*gtd) <= 1e-9
      t = be*t;
    % Safeguard quadratic interpolation
    else
      tq = (-gtd*t^2) / (2*(ft-f-t*gtd));
      if 0.01 <= tq && tq <= 0.99*t 
        t = tq;
      else
        t = be*t;
      end
    end

  end 
  