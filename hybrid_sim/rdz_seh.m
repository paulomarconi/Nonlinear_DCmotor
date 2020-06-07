function [xdot,phi,reset] = rdz(t,x,mode)
% model of the relay switching system with deadzone, d^2x/dt^2 = - rdz(x)
% where rdz(x) = sign(x) for |x| >= 1, 0 for |x| < 1.
%
% "mode" switches when phi = x(1) +/- 1 passes through zero, and that
% instantiates the state event (causes the integrator to switch mode
% from +/-1 to  0 or vice versa).            JH Taylor, 8 March 1996.
%  
%% modified 3-IX-1999 for matlab 5 (changed "if X == []" to "if isempty", 
% added empty output arguments, both to eliminate nasty warning messages)  
%
if isempty(mode), % initialization section
   if abs(x(1)) < 1,
      if x(2) == 0, disp('I will never move ...'), end;
      phi = 0; % => mode = 0
   else
      phi = x(1) - sign(x(1)); % => mode = +/-1
   end
   xdot = []; reset = []; % to prevent new warning msgs!!
   return
end
icall = max(abs(imag(mode)));
if icall == 0, % dynamic equations and switching funs
   if abs(mode) == 1, phi = x(1) - mode;
   elseif mode == 0,
      if abs(x(1)) < 1, phi = 0;
      else phi = x(1) - sign(x(1)); end 
   else disp('mode ~= -1, 0, +1!!'); 
   end
   xdot(1) = x(2);
   xdot(2) = -mode;
   reset = []; % to prevent new warning msgs!!
elseif icall == 1, % corresp to state and/or phi reset
   if abs(real(mode)) == 1, phi = 0; % should => mode = 0 
   else phi = sign(x(1)); end;
   reset = []; xdot = []; % to prevent new warning msgs!!
else
   error('bad value of mode')
end
