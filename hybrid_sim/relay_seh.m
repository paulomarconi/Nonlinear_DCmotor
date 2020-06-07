function [xdot,phi,reset] = relay(t,x,mode)
% model of the relay switching system, d^2x/dt^2 = - sign(x) -- "mode"
% switches when phi = x(1) passes through zero, and that instantiates the
% state event (causes the integrator to switch mode from -1 to 1 or vice
% versa).  JH Taylor, 8 March 1996.
% 
%% modified 3-IX-1999 for matlab 5 (changed "if X == []" to "if isempty",
% added empty output arguments, both to eliminate nasty warning messages) 
%
if isempty(mode), % initialization section
   % Set phi (S) to be consistent with the IC, so that mode
   % will be initialized correctly -- e.g., if x(1) = 0,
   % then x(2) governs mode, i.e. x(2) > 0 => mode = +1 etc.
   if x(1) == 0, phi = x(2);
   else phi = x(1);
   end
   xdot = []; reset = []; % to prevent new warning msgs!!
   return
end
% Define the mode-dependent model and switching flag:
xdot(1) = x(2);
xdot(2) = -mode;
phi = x(1);
reset = []; % to prevent new warning msgs!!
