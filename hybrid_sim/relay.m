function xdot = relay(t,x)
% model of the relay switching system, d^2x/dt^2 = - sign(x) without
% the use of modes and state-event handling.  JH Taylor, 8 March 1996.
% 
xdot(1) = x(2);
xdot(2) = - sign(x(1));
