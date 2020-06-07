function xdot = rdz(t,x)
% model of the relay switching system, d^2x/dt^2 = - rdz(x) where
% rdz(x) = sign(x) for |x| >= 1, 0 for |x| < 1.
%
% Coded to work with ODE45 (no state-event handling; no 'mode')
% JH Taylor, 11 March 1996.
%
if abs(x(1)) > 1, rdz = sign(x(1)); else rdz = 0; end
xdot(1) = x(2);
xdot(2) = -rdz;
