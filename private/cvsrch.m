      function varargout = cvsrch(fun,x,f,g,s,stp,...
                 xtol, maxfev)
%     **********
%
% This is a modified implementation of the algorithm described in:
%
% J. J. More and D. J. Thuente, Line Search Algorithms with Guaranteed
% Sufficient Decrease, TOMS 20 (1994), pp. 286-307.
%
% It attempts to find a step length stp that satisfies the strong Wolfe
% conditions. This function assumes ftol <= gtol and does not used modified 
% function values. 
%
%   Dianne O'Leary   July 1991
%   Yuekai Sun       Mar. 2012
% 
%     **********
      ftol = 0.0001;
      gtol = 0.99;
      stpmin = xtol;
      stpmax = 10;
      
      p5 = 0.5;
      p66 = 0.66;
      xtrapf = 4.;
      info = 0;
      infoc = 1; %#ok<NASGU>
%
%     Compute the initial gradient in the search direction
%
      dginit = g'*s;
%
%     Initialize local variables.
%
      brackt = 0;
      stage1 = 1;
      nfev = 0;
      finit = f;
      dgtest = ftol*dginit;
      width = stpmax - stpmin;
      width1 = 2.*width;
      wa = x;
%
%     The variables stx, fx, dgx contain the values of the step, 
%     function, and directional derivative at the best step.
%     The variables sty, fy, dgy contain the value of the step,
%     function, and derivative at the other endpoint of
%     the interval of uncertainty.
%     The variables stp, f, dg contain the values of the step,
%     function, and derivative at the current step.
%
      stx = 0.0;
      fx = finit;
      dgx = dginit;
      sty = 0.0;
      fy = finit;
      dgy = dginit;
%
%     Start of iteration.
%
   while (1)   
%
%        Set the minimum and maximum steps to correspond
%        to the present interval of uncertainty.
%
         if (brackt) 
            stmin = min(stx,sty);
            stmax = max(stx,sty);
         else
            stmin = stx;
            stmax = stp + xtrapf*(stp - stx);
         end 
%
%        Force the step to be within the bounds stpmax and stpmin.
%
         stp = max(stp,stpmin);
         stp = min(stp,stpmax);
%
%        Evaluate the function and gradient at stp
%        and compute the directional derivative.
%
         x = wa + stp * s;
         if nargout > 6
           [f, g, H] = fun(x);
         else
           [f, g] = fun(x);
         end
         nfev = nfev + 1;
         dg = g' * s;
         ftest1 = finit + stp*dgtest;       
% 
%        Test for convergence.
%
         if (f <= ftest1 && abs(dg) <= gtol*(-dginit)) 
                  info = 1;
         elseif (brackt && stmax-stmin <= xtol) 
                  info = 2;
         elseif (stp == stpmax && f <= ftest1 && dg <= dgtest) 
                  info = 3;
         elseif (stp == stpmin && (f > ftest1 || dg >= dgtest)) 
                  info = 4;
         elseif (nfev >= maxfev) 
                  info = 5;
         end
%
%        Check for termination.
%
         if (info ~= 0) 
                  break
         end
%
%        In the first stage we seek a step for which the modified
%        function has a nonpositive value and nonnegative derivative.
%
         if (stage1 && f <= ftest1 && dg >= min(ftol,gtol)*dginit) 
                stage1 = 0;
         end
% 
%           Call cstep to update the interval of uncertainty 
%           and to compute the new step.
%
         [stx,fx,dgx,sty,fy,dgy,stp,f,dg,brackt,infoc] ...
          = cstep(stx,fx,dgx,sty,fy,dgy,stp,f,dg, ...
                  brackt,stmin,stmax); %#ok<NASGU,ASGLU>
% 
%        Force a sufficient decrease in the size of the
%        interval of uncertainty.
%
         if (brackt) 
            if (abs(sty-stx) >= p66*width1) 
              stp = stx + p5*(sty - stx);
            end
            width1 = width;
            width = abs(sty-stx);
         end
%
%        End of iteration.
%
   end

   if nargout > 6
     varargout = {x,f,g,H,stp,info,nfev};
   else
     varargout = {x,f,g,stp,info,nfev};
   end
%
%     Last card of subroutine cvsrch.
%
      end
