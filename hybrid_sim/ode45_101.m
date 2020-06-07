function [tout, yout] = ode45_101(ypfun, t0, tfinal, y0, tend, tol, trace)
%ODE45_101	Solve differential equations, higher order method.
%	ODE45_101 integrates a system of ODE's using 4th and 5th order 
%	Runge-Kutta formulas, with a hybrid zero-finding algorithm (based on 
%	FZERO) to catch state events. Handles a VECTOR of modes/switching
%	functions, and includes provision for state RESET after each event.
%
%	[T, Y] = ODE45_101('YDOT', T0, Tfinal, Y0) integrates the system
%	of ordinary differential equations described by the M-file YDOT.M,
%	over the interval T0 to Tfinal, with initial conditions Y0.
%	[T, Y] = ODE45('YDOT', T0, Tfinal, Y0, Tend, Tol, Trace) will
%       stop after Tend seconds (useful if your integration is very long
%       and you want to see what's happening), Tol is a non-default value of
%       the tolerance used in the Runge-Kutta-Fehlberg algorithm, and Trace
%       (= 1 or 2)  displays status while the integration proceeds.
%
%       Full calling sequence:
%       [t,y] = ode45_101('ydot',t0,tfinal,y0,tend,tol,trace)
%
%	INPUT:
%	YDOT  - String = name of user-supplied problem description (m-file).
%	        Call: [ydot,phi,reset] = fun(t,y,mode) where:
%	        t      - Time (scalar).
%	        y      - Solution column-vector.
%		mode   - Governs "mode" of dynamic model.
%	        ydot   - derivative column-vector; ydot(i) = dy(i)/dt.
%		phi    - switching function signalling mode change.
%		reset  - reset value of the state y after a state-event.
%	t0    - Initial value of t.
%	tfinal- Final value of t.
%	y0    - Initial value column-vector.
%	tend  - Simulation will stop at tend (sec) (Default: tend = 600).
%	tol   - The desired accuracy. (Default: tol = 1.e-6).
%	trace - If nonzero, each step is printed. (Default: trace = 0).
%
%	OUTPUT:
%	T  - Returned integration time points (column-vector).
%	Y  - Returned solution, one solution column-vector per tout-value.
%
%	The result can be displayed by: plot(tout, yout).
%
%       See also TRAP_101 (lower order/less elegant method).
%
%       D. Kebede, J. H. Taylor,  University of New Brunswick, May, 1995.

%%      Changes by J. H. Taylor, 8 June 1995:
%%      Locked out state-event handling when RK-Fehlberg => point not accepted
%%      Modified logic for iteration for determining time of state event
%%      Modification: modes = 1, 0, -1 not just 1, -1; Initialization changed
%%      also; new conventions are: mode = [] for init, mode = mode + j for reset,
%%      elements of phi must be reset to zero for "on-event" mode (mode = 0).
%%      April 1999 - revised for new matlab syntax requirements; isempty(X) used
%%      to replace "if X == []", and reserved word "catch" replaced with "gotcha".
%%      July 2003 - more syntax changes to accommodate new matlab versions --
%%      changed loop indices from j to jj, other minor nuisance changes...
%%      More changes by JH Taylor, 2 February 2005: 
%% /not yet/ Changed t0, tfinal to tspan, to be consistent with current matlab
%%      Modified initialization logic, to prevent phi(t0) = 0, phi(h) ~=0
%%      from causing an unwanted state event.

% The Fehlberg coefficients:
alpha = [1/4  3/8  12/13  1  1/2]';
beta  = [ [    1      0      0     0      0    0]/4
          [    3      9      0     0      0    0]/32
          [ 1932  -7200   7296     0      0    0]/2197
          [ 8341 -32832  29440  -845      0    0]/4104
          [-6080  41040 -28352  9295  -5643    0]/20520 ]';
gamma = [ [902880  0  3953664  3855735  -1371249  277020]/7618050
          [ -2090  0    22528    21970    -15048  -27360]/752400 ]';
pow = 1/5;
if nargin < 5, tend = 300; end %%% five minutes max!
if nargin < 6, tol = 1.e-6; end
if nargin < 7, trace = 0; end

% Initialization
t = t0;
hmax = (tfinal - t)/16;
h = hmax/8;
y = y0(:);
%%%%%%% chunk = 256;
%%%%%%% tout = zeros(chunk,1);
%%%%%%% yout = zeros(chunk,length(y));
k = 1;
tout(k) = t;
yout(k,:) = y.';
j_imag = j; 	% initialize "j" for reset

% Call model so that it can provide phi for initial mode evaluation:
mode = [];
[junk,phi,trash] = feval(ypfun,t,y,mode);
mode = sign(phi);
mdim = length(mode);
% initialize arrays for handling RKF points and possible multiple state events
f = zeros(length(y),6);
ym = zeros(length(y),mdim + 1);
tv = zeros(1,mdim + 1);
phim = zeros(mdim, mdim + 1);
cv = zeros(1,mdim + 1);
% set the elapsed time stopwatch:
timer = clock;
% also, set the total time stopwatch:
max_timer = clock;
if trace
   msg = ['Point A (init):  t_init = ',num2str(t),' ,  h = ',num2str(h)];
   clc, disp(msg)
   if trace == 2, y_0 = y', mode_0 = mode, end
end

% The main loop
while (t < tfinal) & (t + h > t)
   if t + h > tfinal, h = tfinal - t; end
   % report progress in simulation (every 30 seconds): 
   if etime(clock,timer) > 30, 
      timer = clock; 
      tot_elapsed = etime(clock,max_timer);
      disp(['T_real = ',num2str(tot_elapsed),' ,  T_sim = ',num2str(t)]);
   end 
   % see if it is time to give up completely:
   if etime(clock,max_timer) > tend, 
      simu_time = t, break
   end 
   % Save the last "accepted" point in case of state event:
   yold = y(:);  told = t;   phiold = phi(:);
   % At the start of each step the new point is NOT accepted until:
   % (1) the RKF error handler is satisfied, (2) the SEH is satisfied
   accepted = 0;
   % Compute the derivatives
   y = y(:);
   [temp,junk1,junk2] = feval(ypfun,t,y,mode);
   f(:,1) = temp(:);
   fold(:,1) = f(:,1);
   for jj = 1:5
      [temp,junk1,junk2] = feval(ypfun, t+alpha(jj)*h, y+h*f*beta(:,jj), mode);
      f(:,jj+1) = temp(:); 
   end
   % Estimate the error and the acceptable error
   delta = norm(h*f*gamma(:,2),'inf');
   tau = tol*max(norm(y,'inf'),1.0);
   % Update the solution only if the error is acceptable
   if delta <= tau
      accepted = 1;
      t = t + h;
      y = y + h*f*gamma(:,1); %%% y = acceptable y_RKF solution
      % Prepare for SE test: Get the current trial point phi:
      [junk1,phi,junk2] = feval(ypfun, t, y, mode); phi=phi(:);
      % load SEH arrays:
      for im = 1:mdim+1
         ym(:,im) = y;
         tv(im) = t;
	 phim(:,im) = phi;
         cv(im) = 0; 
      end
      if trace
         msg = ['Point B (y_RKF):  t = ',num2str(t),' ,  h = ',num2str(h)];
         disp(msg)
         if trace == 2, y_RKF = y', mode_RKF = mode, end
      end
   end
   % Update the step size
   if delta ~= 0.0
      h = min(hmax, 0.8*h*(tau/delta)^pow);
   end

   gotcha = 0;
   % Now: we only test and handle state events when RKF => acceptable:
   if ( (accepted == 1) & (k > 1) ),
      for im =1:mdim 
	if sign(phi(im)) ~= sign(phiold(im)) % NEW TEST
          %% OK - we have detected a state event:
	  cv(im) = 1;  %% point IS a state-event
          niter = 0;
	  fini = 1;
	  % NB: bias phi so that the boundary is crossed for sure (biasing
	  % phi using phi(im) works for either type of state event):
	  bias = -eps*sign(phi(im));
	  ha=0; hb=h; hc=0; hd=h; he=h; % init hd & he to elim complaints...
	  phia = phiold(im)+bias; phib = phi(im)+bias; phic = phiold(im)+bias;
	  % Main loop, exit from middle of the loop
	  while (fini > 0)
            % hybrid zero-finding algorithm based on FZERO.M
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
            hm = 0.5*(hc - hb);
            % Convergence test and possible exit based on "hb"
            toler = 2*eps*max(abs(hb),1.0);
            if (abs(hm) <= toler) | (phib == 0)
               % NB: abs(hm) <= eps => as close as you're gonna get!
               % see which "h" to use to make sure SE has happened:
	       if sign(phib) ~= sign(phiold(im) + bias), hb = hb;
	         else, hb = hc;
	       end
	       f(:,1) = fold(:,1);
	       for jj = 1:5
		  [temp,junk1,junk2] = feval(ypfun, told+alpha(jj)*hb, ...
                    yold+hb*f*beta(:,jj), mode);
		  f(:,jj+1) = temp(:);
	       end
	       tstar = told + hb;
	       ystar = yold + hb*f*gamma(:,1); % y(+) or y(-) depending on exit
	       [junk,phi_t,trash] = feval(ypfun, tstar, ystar, mode); phi_t = phi_t(:); 
	       tv(im) = tstar;
	       ym(:,im) = ystar;
	       phim(:,im) = phi_t;
	       break,
            end
            % Choose bisection or interpolation
            if (abs(he) < eps) + (abs(phia) <= abs(phib))
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
               if (2.0*zp < 3.0*hm*zq - abs(eps*zq)) * (zp < abs(0.5*he*zq))
                  he = hd;  hd = zp/zq;
               else
                  hd = hm;  he = hm;
               end;
            end % bisection/interpolation
            % evaluate next point
            ha = hb;
            phia = phib;
            if abs(hd) > eps, hb = hb + hd;
            else if hb > hc, hb = hb - eps;
                 else hb = hb + eps;
                 end
            end
	    hstar = hb;   % use new "h" to compute the slopes again
	    f(:,1) = fold(:,1);
	    for jj = 1:5
		[temp,junk1,junk2] = feval(ypfun, told+alpha(jj)*hstar, ...
                    yold+hstar*f*beta(:,jj), mode);
		f(:,jj+1) = temp(:);
	    end
	    tstar = told + hstar; 
	    ystar = yold + hstar*f*gamma(:,1); % y(+) or y(-) depending on exit
	    [junk,phi_t,trash] = feval(ypfun, tstar, ystar, mode);
	    phi_t = phi_t(:);  phib = phi_t(im) + bias;
            % normal exit (converged to solution):
            if (abs(phib) < eps),
		tv(im) = tstar; 
		ym(:,im) = ystar;
		phim(:,im) = phi_t;
		fini = 0;
            end
	  end  % of WHILE (fini > 0) loop
          if trace 
            msg = [ 'N_iter = ',num2str(niter),' im = ',num2str(im), ...
                  ' Catch = ',num2str(cv(im)),' t_se = ',num2str(tv(im)) ];
            disp(msg);
            time1 = tstar, y_rk1= ystar'
          end
	end  % of IF sign(phi(im)) ~= sign(phiold(im)) statement
      end  % of FOR im = 1:mdim loop
      % Now, check all modes for possible/earliest change:
      tmin = t;
      imin = mdim + 1;   %% values from RKF integration
      for im = 1:mdim
         if (tv(im) < tmin) | ((abs(tv(im) - tmin) < eps) & cv(im) == 1),
            tmin = tv(im);
            imin = im;
	 else
	    cv(im) = 0; % discard later ones (preliminary)
         end
      end
      t = tmin;   % this is a bone fide point
      y = ym(:,imin);
      phi = phim(:,imin);
      gotcha = cv(imin);
      if gotcha == 1,
	 % one or more SEs have happened:
         if trace
            msg=['Point C (SEH):  t = ',num2str(t), ...
               ' ,  imin = ',num2str(imin)];
            disp(msg)
            if trace == 2, y_SEH = y', mode_SEH = mode, phi_seh = phi, end
         end
	 % first, check for simultaneous SEs:
	 num_ses = 0; 
	 for im = 1:mdim
	    if sign(phim(im,im)) == sign(phim(im,imin)) & cv(im) == 1,
	     num_ses = num_ses + 1;
	    else
	     cv(im) = 0;
	    end
	 end
	 if num_ses == 0, disp(['Oops: SE glitch at t = ',num2str(tv(im))]); end
	 if num_ses > 1, disp(['Note: simultaneous SEs at t = ',...
             num2str(tv(im)), ' Number of SE = ',num2str(num_ses)]); end
	 % get the reset value of the state (make abs(mode(imin)) = 2):
	 % (NB: make mode(im s.t. cv(im) = 1) complex):
	 rmode = mode;
	 for im = 1:mdim
	    if cv(im) == 1, rmode(im) = rmode(im) + j_imag; end
	 end
	 % (we assume the model can handle multiple simult. SEs!)
	 [junk,phi_reset,y_reset] = feval(ypfun, t, y, rmode);
	 % see if reset is called for:
	 if ~isempty(y_reset),
	     % put the captured point in the output (make room for state reset)
	     k = k+1;
%%%%%%% 	     if k > length(tout)
%%%%%%% 		tout = [tout; zeros(chunk,1)];
%%%%%%% 		yout = [yout; zeros(chunk,length(y))];
%%%%%%% 	     end   
	     tout(k) = t;
	     yout(k,:) = y.';
	     % set up reset value to be put on next:
	     y = y_reset(:);   %%% y_reset
	     if trace
		msg = ['Point D (reset):  t = ',num2str(t)];
		disp(msg)
		if trace == 2, y_reset = y', mode_reset = mode, end
	     end
	 end % of y_reset
	 % now, at last, update mode:
	 for im = 1:mdim
	   if cv(im) == 0, % nothing to be done
	   % see if we are entering "on event" situation:
	   elseif phi_reset(im) == 0, phi(im) = 0; mode(im) = 0;
	   % otherwise implement crossing event:
	   else
	      if sign(phi(im)) == sign(mode(im)),
		 t, y_rkf = y', phi_rkf = phi', mode, im,
		 disp('mode and phi disagree!'); 
	      end
	      mode(im) = sign(phi(im));
	   end
	 end % of mode updating
      end % of IF gotcha == 1 statement
   end  % of IF accepted == 1 statements
      % OK: we must have a good point to add - y_RKF, y_se, or y_reset:
      k = k+1;
%%%%%%%       if k > length(tout)
%%%%%%%          tout = [tout; zeros(chunk,1)];
%%%%%%%          yout = [yout; zeros(chunk,length(y))];
%%%%%%%       end
      tout(k) = t;
      yout(k,:) = y.';
      if trace, disp('Point E:  End of state-event handling'); end
end  % of WHILE (t < tfinal) & (t + h > t) loop
if trace, disp('Point F:  End of integration loop'); end

if (t < tfinal)
   disp(['Singularity likely, or time ran out (tend = ',num2str(tend),')'])
   t
end

tout = tout(1:k);
yout = yout(1:k,:);
