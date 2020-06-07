function [xdot,phi,reset] = twin_ball(t,x,mode)
%
% design convention: mode = [] for init, = -1, 0, 1 for different dynamic
% equations, = -1+j, 1+j for state reset prior to mode = 1, 0, -1.  Note:
% if ANY element of mode is complex then reset it (only).  JH Taylor 17-X-1996
%
%% modified 3-IX-1999 for matlab 5 (changed "if X == []" to "if isempty",
% added empty output arguments, both to eliminate nasty warning messages) 
%
% at initialization mode = []; return -> length and value of mode vector.
if isempty(mode),
   for i=1:2,
      if x(2*i-1) == 0, phi(i) = x(2*i);
      else phi(i) = x(2*i-1);
      end
   end
   xdot = []; reset = []; % to prevent new warning msgs!!
   return
end
%
icall = max(abs(imag(mode)));
if icall == 0, % mode real -> dynamic equations and switching funs
   xdot(1) = x(2);
   xdot(2) = -mode(1);
   xdot(3) = x(4);
   xdot(4) = -mode(2);
   phi(1) = 100*x(1);
   phi(2) = 100*x(3);
   reset = []; % to prevent new warning msgs!!
elseif icall == 1, % mode complex corresp to state reset
   phi = -mode;
   reset = x; 
   for i=1:2,  % only reset states for the actual mode change:
      if abs(imag(mode(i))) == 1,
         reset(2*i-1) = 0;
         reset(2*i) = 0.8*x(2*i);
      end
   end
   xdot = []; % to prevent new warning msgs!!
else
   error('bad value of mode')
end
