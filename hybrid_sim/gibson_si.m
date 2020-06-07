function [xdot,phi,reset] = relay(t,x,mode)
%% version of Gibson roll-control problem suitable for
%% simulation via ode45_101 (using modes)
%
% JH Taylor - 24 March 1999

%%% converted to SI units
% states x1 = phi, x2 = phidot, x3, x4 = servoamp model states
% define the servoamp dynamics (convert transfer fn to state space)
% A = [ 0 1 ; -2000 -120 ]; B = [ 0 ; 1 ]; C = [ 2000 0 ]; D = 0;

h = 22.24; F0 = 444.8; % thruster "relay" constants (Newtons)
J = 4.68; % N-m/sec^2; = 3.45 lb-ft/sec^2
R0 = 0.6096; % m; = 2 ft
Kp = 1868; % N/rad; = 420 lb/rad
Kv = 186.8; % N/(rad/sec); = 42 lb/(rad/sec)

% initializes the mode in ode45_101
if isempty(mode);
  e = 2000*x(3);
  if e > h, phi(1) = 1.0;
  elseif e < -h, phi(1) = -1.0;
  else error('sorry - cannot initialize properly . . .');
  end;
  xdot = [];
  reset = [];
  return;
end;

% actual model and switching logic
e = 2000*x(3);
icall = max(abs(imag(mode)));
if icall == 0, % dynamic equations and switching funs
   xdot(1) = x(2);
   thrust = F0 * mode;
   xdot(2) = R0*thrust/J;
   xdot(3) = x(4);
   Uamp = - Kv*x(2) - Kp*x(1);
   xdot(4) = -2000*x(3) -120*x(4) + Uamp;
   if mode == 1, phi(1) = e + h;
   elseif mode == -1, phi(1) = e - h;
   else error('sorry - mode should never be 0!!');
   end
   reset = [];
elseif icall == 1, % corresp to possible state reset
   xdot = [];
   phi(1) = NaN;
   reset = []; 
   disp(['switching at t = ',num2str(t)])
else
   error('bad value of mode')
end
% end of relay with hysteresis model
