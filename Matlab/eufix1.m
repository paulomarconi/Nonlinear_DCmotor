function [tout, xout] = eufix1(dxfun, tspan, x0, stp, trace)
%EUFIX1	Solve ordinary state-vector differential equations, low order method.
%	EUFIX1 integrates a set of ODEs xdot = f(x,t) using the most
%	elementary Euler algorithm, without step-size control.
%
%	CALL:
%       [t, x] = eufix1('dxfun', tspan, x0, stp, trace)
%
%	INPUT:
%	dxfun - String containing name of user-supplied problem description.
%	        Call: xdot = model(t,x) coded in fname.m => dxfun = 'fname'.
%	        t     - Time (scalar).
%	        x     - Solution column-vector at time t.
%	        xdot  - Returned derivative column-vector; xdot = dx/dt.
%	tspan - Range of t for the desired solution; tspan = [t0 tf].
%	tf    - Final value of t.
%	x0    - Initial value column-vector.
%   	stp  - The specified integration step (default: stp = 1.e-2).
%	trace - If nonzero, each step is printed (default: trace = 0).
%
%	OUTPUT:
%	t  - Returned integration time points (row-vector).
%	x  - Returned solution, one column-vector per tout-value.
%
%	Display result by: plot(t, x) or plot(t, x(:,2)) or plot(t, x(:,2), x(:,5)).

% Initialization
if nargin < 4, stp = 1.e-2; disp('H = 0.02 by default'); end
if nargin < 5, trace = 0; end     %% disable trace if not requested
t0 = tspan(1); tf = tspan(2);
if tf < t0, error('tf < t0!'); return; end  %% check for glaring error
t = t0;
h = stp;
x = x0(:);
k = 1;
tout(k) = t;      % initialize output arrays
xout(k,:) = x.';
if trace
   clc, t, h, x
end

% The main loop

while (t < tf) 
   if t + h > tf, h = tf - t; end
   % Compute the derivative
   dx = feval(dxfun, t, x); dx = dx(:);
   % Update the solution (with no check on error)
   t = t + h;
   x = x + h*dx;
   k = k+1;
   tout(k) = t;
   xout(k,:) = x.';
   if trace
      home, t, h, x, dx
   end
end
if (t < tf) % if true, something bad happened!
   disp('Singularity or modeling error likely.')
   t
end
% ... here is the output (tout in row vector form)
tout = tout(1:k);
xout = xout(1:k,:);
