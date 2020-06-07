function [tout, yout] = trap_101(dyfun, t0, tf, y0, step, tol, trace)
%  TRAP_101   Solve differential equations, "trapezoidal" method.
%       mod 6 -- 4 March 1996
%	TRAP_101 integrates a system of ordinary differential equations
%	using the "trapezoidal" algorithm, with a hybrid interpolation
%       scheme to catch state events (points where ydot is discontinuous).
%
%       CALL: [tout,yout] = trap_101(dyfun,t0,tf,y0,step,tol,trace)
%	INPUT:
%	dyfun - String containing name of user-supplied problem description.
%	        Call: [ydot,phi,reset] = fun(t,y,mode) where dyfun = 'fun'.
%	        t     - Time (scalar).
%	        y     - Solution column-vector.
%		mode  - Governs "mode" of the dynamic model
%	        ydot  - Returned derivative column-vector; ydot(i) = dy(i)/dt.
%		phi   - Returned switching function signalling mode change
%		reset - Returned reset value of the state y
%	t0    - Initial value of t.
%	tf    - Final value of t.
%	y0    - Initial value column-vector.
%   	step  - The specified integration step. (Default: step = 1.e-2).
%	tol   - Zero-finding tolerance. (Default: tol = eps).
%	trace - If nonzero, each step is printed, if ==2 zero-finding debug
%		output is provided, if == 3 state-event handling debug out-
%		put is provided, if == 4 you get both. (Default: trace = 0).
%
%	OUTPUT:
%	T  - Returned integration time points (column-vector).
%	Y  - Returned solution, one solution column-vector per state.
%
%% JH Taylor, University of New Brunswick, 13 May 1995 and subsequent mods.
%% April 1999 - revised for new matlab syntax requirements; isempty(X) used  
%% to replace "if X == []", and reserved work "catch" replaced with "gotcha".

% This routine contains part of MATLAB's FZERO algorithm, with permission.
%% Initialization
if nargin < 5, step = 1.e-2; end
if nargin < 6, tol = eps; end
if nargin < 7, trace = 0; end
disp('This is routine TRAP_101, version 6; 09 April 1999');
t = t0;
h = step;
y = y0(:);
%% Call model to get phi for initial mode evaluation:
mode = [];
[junk,phi,trash] = feval(dyfun,t,y,mode);
mode = sign(phi);
mdim = length(mode);
%% allocate arrays, store initial condition:
chunk = 256;
tout = zeros(chunk,1);
yout = zeros(chunk,length(y));
k = 1;
tout(k) = t;
yout(k,:) = y.';
% initialize arrays for handling state events
ym = zeros(length(y),mdim + 1);
phim = zeros(mdim,mdim + 1);
tv = zeros(1,mdim + 1);
cv = zeros(1,mdim + 1);
% set the elapsed time stopwatch:
timer = clock;
% also, set the total time stopwatch:
max_timer = clock;
if trace, dt = h, time_0 = t, y_0 = y', phi_0 = phi', m_0 = mode', end
%
%% Main integration loop
while (t < tf) & (t + h > t)
   if t + h > tf, h = tf - t; end
   % report progress in simulation:
   if (trace & etime(clock,timer) > 30),
      timer = clock;
      tot_elapsed = etime(clock,max_timer);
      disp(['T_real = ',num2str(tot_elapsed),' ,  T_sim = ',num2str(t,6)]);
   end
   % save data from last "accepted" point:
   yold = y(:); told = t; phiold = phi;
   % Also: compute the slope at the accepted point:
   [ydotold,junk,trash]=feval(dyfun,t,y,mode); ydotold=ydotold(:);
   % Update the solution (TRIAL point!)
   tt = t + h;
   yp = yold + h*ydotold;  % use old ydot to predict; get new ydot ...
   [ydot,junk,trash] = feval(dyfun,tt,yp,mode); ydot = ydot(:);
   yt = yold + h*(ydotold + ydot)/2; % ... and correct (trapezoidal rule)
   % Update the (trial point) phi:
   [junk,phi,trash] = feval(dyfun,tt,yt,mode); phi=phi(:);
   % load into arrays for use in SEH:
   for im = 1:mdim+1
      ym(:,im) = yt;
      phim(:,im) = phi;
      tv(im) = tt;
      cv(im) = 0;
   end
   gotcha = 0;
   for im = 1:mdim
      %% if phi(im)*phiold(im) < 0, THIS TEST IS BAD FOR MODE == 0 CASES
      if sign(phi(im)) ~= sign(phiold(im)) % NEW TEST
         %% possible state event detected:
         gotcha = 1;
         cv(im) = 1;
         niter = 0;
         % set up fzero-like zero-finding scheme:
         fini = 1;
         % NB: bias phi so that the boundary is surely crossed:
         bias = -tol*sign(phi(im));
         phia = phiold(im)+bias; phib = phi(im)+bias; phic = phia;
         ha=0; hb=h; hc=0; hd=h; he=h; % init hd & he to elim complaints...
         while fini > 0
            % hybrid zero-finding algorithm based on fzero
            % hb = best h so far, ha is previous hb, and hc is an h
            % such that phi has the opposite sign from that at hb.
            niter = niter + 1;
            if (phib > 0) == (phic > 0)
               hc = ha;  phic = phia;
               hd = hb - ha;  he = hd;
            end
            if abs(phic) < abs(phib)
               ha = hb;  hb = hc;  hc = ha;
               phia = phib;  phib = phic;  phic = phia;
            end
            if (trace  == 2) | (trace  == 4),
               msg1 = ['ha  =  ',num2str(ha,8),' hb  =  ',num2str(hb,8),...
                  ' hc  =  ',num2str(hc,8),' hd  =  ',num2str(hd,8),...
                  ' he  =  ',num2str(he,8)];
               msg2 = ['phia = ',num2str(phia,8),' phib = ',num2str(phib,8),...
                  ' phic = ',num2str(phic,8)];
               disp(msg1);
               disp(msg2);
            end
            hm = 0.5*(hc - hb);
            % Convergence test and possible side-door exit:
            toler = 2.0*tol*max(abs(hb),1.0);
            if (abs(hm) <= toler | phib == 0.0)
               % NB: abs(hm) <= toler => as close as you're gonna get!
               % see which "h" to use to make sure SE has happened:
               if sign(phib) ~= sign(phiold(im) + bias), hb = hb;
                  else hb = hc; end % hb or hc _has_ to be past the SE
               %% evaluate the current y and phi:
               tb = told + hb;
               yp = yold + hb*ydotold;
               [ydot,junk,trash] = feval(dyfun,tb,yp,mode); ydot = ydot(:);
               yb = yold + hb*(ydotold + ydot)/2; % trapezoidal rule again
               [junk,phit,trash] = feval(dyfun,tb,yb,mode); phit=phit(:); 
               tv(im) = tb;
               ym(:,im) = yb;
               phim(:,im) = phit;
               fini = 0;
            end
            % Choose bisection or interpolation
            if (abs(he) < tol) + (abs(phia) <= abs(phib))
               % Bisection
               hd = hm;  he = hm;
            else
               % Interpolation
               zs = phib/phia;
               if (ha == hc)
                  % Linear interpolation
                  zp = 2.0*hm*zs;
                  zq = 1.0 - zs;
               else
                  % Inverse quadratic interpolation
                  zq = phia/phic;
                  zr = phib/phic;
                  zp = zs*(2.0*hm*zq*(zq - zr) - (hb - ha)*(zr - 1.0));
                  zq = (zq - 1.0)*(zr - 1.0)*(zs - 1.0);
               end;
               if zp > 0, zq = -zq; else zp = -zp; end;
               % Is interpolated point acceptable
               if (2.0*zp < 3.0*hm*zq - abs(tol*zq)) * (zp < abs(0.5*he*zq))
                  he = hd;  hd = zp/zq;
               else
                  hd = hm;  he = hm;
               end;
            end % bisection/interpolation
            % evaluate next point
            ha = hb;
            phia = phib;
            if abs(hd) > tol, hb = hb + hd;
            else if hb > hc, hb = hb - tol;
                 else hb = hb + tol;
                 end
            end
            tb = told + hb;
            yp = yold + hb*ydotold;
            [ydot,junk,trash] = feval(dyfun,tb,yp,mode); ydot = ydot(:);
            yb = yold + hb*(ydotold + ydot)/2; % trapezoidal rule again
            [junk,phit,trash] = feval(dyfun,tb,yb,mode); phit = phit(:);
            phib=phit(im) + bias;
            % normal exit (converged to solution):
            if abs(phib) < tol
               tv(im) = tb;
               ym(:,im) = yb;
               phim(:,im) = phit;
               fini = 0;
            end
         end % of WHILE fini > 0 loop
         if (trace  == 2) | (trace  == 4),
            msg = ['N_iter = ',num2str(niter),' im = ',num2str(im),' cv(im) = ',...
               num2str(cv(im)),' tv(im) = ',num2str(tv(im)) ];
            disp(msg) 
            time_seh = tb, y_seh = yb'
         end
      end % of IF sign(phi(im)) ~= sign(phiold(im)) statement
   end % of FOR im = 1:mdim loop
   if gotcha == 0 % we didn't catch a SE
      imin = mdim + 1;  % the trial point is accepted
      t = tv(imin);
      y = ym(:,imin);
      phi = phim(:,imin);
   else
      % check for earliest change:
      tmin = tt; %% values from trapezoidal integration
      imin = mdim + 1;
      for im = 1:mdim
         if (tv(im) < tmin), 
            tmin = tv(im);
            imin = im;
         end
      end
      t = tmin;  % this is an accepted point (with SE(s))
      y = ym(:,imin);
      phi = phim(:,imin);
      if (trace == 3) | (trace  == 4),
         msg=['Point C (SEH):  t = ',num2str(t), ...
            ' ,  imin = ',num2str(imin)];
         disp(msg)
         y_SEH = y', phi_SEH = phi', cv_SEH = cv, mode_SEH = mode
      end
      % one or more SEs have happened; rigorous test for simultaneity:
      num_ses = 0;
      for im = 1:mdim
         if (sign(phim(im,imin)) ~= sign(phiold(im)) & cv(im) == 1),
            num_ses = num_ses + 1;
            %% cv(im) = 1; % be sure to retain it
         else
            cv(im) = 0; % discard it
         end
      end
      if (trace & num_ses == 0),
         disp(['Oops: SE glitch at t = ',num2str(tv(im))]); end
      if (trace & num_ses > 1), disp(['Note: simultaneous SEs at t = ',...
         num2str(tv(im))]); end
      % get the reset value of the state and/or of phi
      % (NB: make mode(im s.t. cv(im) = 1) complex):
      rmode = mode;
      for im = 1:mdim
         if cv(im) == 1, rmode(im) = rmode(im) + j; end
      end
      [junk,phi_reset,y_reset] = feval(dyfun, t, y, rmode);
      % see if state reset is called for:
      if ~isempty(y_reset),
         % put the captured point in the output (make room for state reset)
         k = k+1;
         if k > length(tout)
            tout = [tout; zeros(chunk,1)];
            yout = [yout; zeros(chunk,length(y))];
         end   
         tout(k) = t;
         yout(k,:) = y.';
         % set up state reset value to be put on next:
         y = y_reset(:); %% y_reset; time doesn't change
         if (trace == 3) | (trace  == 4),
            msg = ['Point D (y_reset):  t = ',num2str(t)];
            disp(msg)
            y_reset = y', phi_reset = phi', cv_reset = cv, mode_reset = mode
         end   
      end % of y_reset
      % now, at last, update the mode:
      for im = 1:mdim
         if cv(im) == 0, % nothing to be done
         % see if we are entering "on-boundary" situation:
         elseif phi_reset(im) == 0, phi(im) = 0; mode(im) = 0;
         % otherwise implement crossing event:
         else 
            if sign(phi(im)) == sign(mode(im)),
               disp(['Inconsistent mode at t = ',num2str(tv(im))]); end
            mode(im) = sign(phi(im));
         end
      end % of mode updating
      if (trace == 3 | trace == 4), mode_update = mode, end
   end % of IF catch == 0 statement
   % we have a good point to add - y_trap, y_se, or y_reset:
   k = k+1;
   if k > length(tout)
      tout = [tout; zeros(chunk,1)];
      yout = [yout; zeros(chunk,length(y))];
   end   
   tout(k) = t;
   yout(k,:) = y.';
   if trace, time_k = t, y_k = y', end
end % of main WHILE loop
if (t < tf)
   disp(['Premature end of simulation, t = ',num2str(t)]);
end
tout = tout(1:k);
yout = yout(1:k,:);
