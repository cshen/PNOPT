     function  [stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,info] ...
       = MoreThuenteStep(stx,fx,dx,sty,fy,dy,stp,fp,dp,brackt,stpmin,stpmax)
%   **********
%
% MoreThuenteStep seeks a safeguarded step for a line search and updates an 
% interval of uncertainty that brackets a minimizer of the function. 
%
%   Dianne O'Leary   July 1991
%   Yuekai Sun       Mar. 2012
% 
%   **********
      p66 = 0.66;
      info = 0;
%
%     Check the input parameters for errors.
%
%       if ((brackt && (stp <= min(stx,sty) || ...
%           stp >= max(stx,sty))) || ...
%           dx*(stp-stx) >= 0.0 || stpmax < stpmin) 
%          return
%       end
%
%     Determine if the derivatives have opposite sign.
%
      sgnd = dp*(dx/abs(dx));
      
      if (isinf(fp) || isnan(fp)) || (isinf(dp) || isnan(dp))
        bound = 0;
        stpf = 0.5*(stx+sty);
%
%     First case. A higher function value.
%     The minimum is bracketed. If the cubic step is closer
%     to stx than the quadratic step, the cubic step is taken,
%     else the average of the cubic and quadratic steps is taken.
%     
      elseif (fp > fx) 
         info = 1;
         bound = 1;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
         gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
         if (stp < stx) 
             gamma = -gamma;
         end
         p = (gamma - dx) + theta;
         q = ((gamma - dx) + gamma) + dp;
         r = p/q;
         stpc = stx + r*(stp - stx);
         stpq = stx + ((dx/((fx-fp)/(stp-stx)+dx))/2)*(stp - stx);
         if (abs(stpc-stx) < abs(stpq-stx)) 
            stpf = stpc;
         else
           stpf = stpc + (stpq - stpc)/2;
         end 
%          brackt = 1;
%
%     Second case. A lower function value and derivatives of
%     opposite sign. The minimum is bracketed. If the cubic
%     step is closer to stx than the quadratic (secant) step, 
%     the cubic step is taken, else the quadratic step is taken.
%
      elseif (sgnd < 0.0) 
         info = 2;
         bound = 0;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
         gamma = s*sqrt((theta/s)^2 - (dx/s)*(dp/s));
         if (stp > stx) 
            gamma = -gamma;
         end
         p = (gamma - dp) + theta;
         q = ((gamma - dp) + gamma) + dx;
         r = p/q;
         stpc = stp + r*(stx - stp);
         stpq = stp + (dp/(dp-dx))*(stx - stp);
         if (abs(stpc-stp) > abs(stpq-stp))
            stpf = stpc;
         else
            stpf = stpq;
         end 
%          brackt = 1;
%
%     Third case. A lower function value, derivatives of the
%     same sign, and the magnitude of the derivative decreases.
%     The cubic step is only used if the cubic tends to infinity 
%     in the direction of the step or if the minimum of the cubic
%     is beyond stp. Otherwise the cubic step is defined to be 
%     either stpmin or stpmax. The quadratic (secant) step is also 
%     computed and if the minimum is bracketed then the the step 
%     closest to stx is taken, else the step farthest away is taken.
%
      elseif (abs(dp) < abs(dx)) 
         info = 3;
         bound = 1;
         theta = 3*(fx - fp)/(stp - stx) + dx + dp;
         s = norm([theta,dx,dp],inf);
%
%        The case gamma = 0 only arises if the cubic does not tend
%        to infinity in the direction of the step.
%
         gamma = s*sqrt(max(0.,(theta/s)^2 - (dx/s)*(dp/s)));
         if (stp > stx) 
             gamma = -gamma;
         end
         p = (gamma - dp) + theta;
         q = (gamma + (dx - dp)) + gamma;
         r = p/q;
         if (r < 0.0 && gamma ~= 0.0)
            stpc = stp + r*(stx - stp);
         elseif (stp > stx)
            stpc = stpmax;
         else
            stpc = stpmin;
         end 
         stpq = stp + (dp/(dp-dx))*(stx - stp);
         if (brackt) 
            if (abs(stp-stpc) < abs(stp-stpq)) 
               stpf = stpc;
            else
               stpf = stpq;
            end
         else
            if (abs(stp-stpc) > abs(stp-stpq)) 
               stpf = stpc;
            else
               stpf = stpq;
            end 
         end 
%
%     Fourth case. A lower function value, derivatives of the
%     same sign, and the magnitude of the derivative does
%     not decrease. If the minimum is not bracketed, the step
%     is either stpmin or stpmax, else the cubic step is taken.
%
      else
         info = 4;
         bound = 0;
         if (brackt) 
            theta = 3*(fp - fy)/(sty - stp) + dy + dp;
            s = norm([theta,dy,dp],inf);
            gamma = s*sqrt((theta/s)^2 - (dy/s)*(dp/s));
            if (stp > sty) 
                gamma = -gamma;
            end
            p = (gamma - dp) + theta;
            q = ((gamma - dp) + gamma) + dy;
            r = p/q;
            stpc = stp + r*(sty - stp);
            stpf = stpc;
         elseif (stp > stx)
            stpf = stpmax;
         else
            stpf = stpmin;
         end 
      end 
%
%     Update the interval of uncertainty. This update does not
%     depend on the new step or the case analysis above.
%
      if (fp > fx) 
         sty = stp;
         fy = fp;
         dy = dp;
         brackt = 1;
      else
         if (sgnd < 0.0)
            sty = stx;
            fy = fx;
            dy = dx;
            brackt = 1;
         end 
         stx = stp;
         fx = fp;
         dx = dp;
      end
%
%     Compute the new step and safeguard it.
%
      stpf = min(stpmax,stpf);
      stpf = max(stpmin,stpf);
      stp = stpf;
      if (brackt && bound)
         if (sty > stx) 
            stp = min(stx+p66*(sty-stx),stp);
         else
            stp = max(stx+p66*(sty-stx),stp);
         end
      end
      return
%
%     Last card of subroutine cstep.
%
     end