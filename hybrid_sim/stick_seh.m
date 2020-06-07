function [xdot,phi,reset] = stiction(t,x,mode)
% model of an electro-mechanical system with stiction and saturation
%
% "mode" switches from +/- 1 to 0 when velocity passes through zero, and 
% transitions back to +/- 1 when there is sufficient torque to overcome
% stiction.  States are x(1) = position, x(2) = velocity; the model is
% identical to Taylor & Strobel, 1985 ACC paper.  JH Taylor, 8-VI-1996
%% modified 3-IX-1999 for matlab 5 (changed "if X == []" to "if isempty",
% added empty output arguments, both to eliminate nasty warning messages)
%
global v_in omega
if isempty(mode), % initialization section
   if abs(x(2)) == 0,
      disp('Initialized stuck ...');
      phi = 0; % => mode = 0
   else
      disp('Initialized moving ...');
      phi = sign(x(2)); % => mode = +/-1
   end
   xdot = []; reset = []; % to prevent new warning msgs!!
   return
end
m_1 = 5.0;    % Nm/v
delta = 0.5;  % v 
m_2 = 1.0;    % Nm/v
f_v = 0.1;    % Nms/rad
f_c = 1.0;    % Nm
MoI = 0.01;   % kg-m^2
% first, prelim calcs: 
% saturate the input voltage -> electrical torque:
volts = v_in*sin(omega*t);
av = abs(volts); sv = sign(volts);
if av < delta, T_e = m_1*volts;
else T_e = (m_1*delta + m_2*(av-delta))*sv; end 
% now, either execute the model or reset function:
icall = max(abs(imag(mode)));
if icall == 0, % models for dynamic equations and switching funs
   if abs(mode) == 1,  % motor moving:
       xdot(1) = x(2);
       xdot(2) = (T_e - f_v*x(2) - f_c*mode)/MoI;
       phi = x(2);
   else,               % motor stuck:
       xdot(1) = 0.0;
       xdot(2) = 0.0;
       % NB! phi == 0 until torque suffices for break-away
       if (abs(T_e) - f_c) < 0, phi = 0;
       else phi = (abs(T_e) - f_c)*sign(T_e); end
   end
   reset = []; % to prevent new warning msgs!!
elseif icall == 1, % corresp to state and/or phi reset
   if abs(real(mode)) == 1
       % we are going into "stuck" mode...
       reset(1) = x(1);
       reset(2) = 0.0;  % do this so we EXACTLY stick
       phi = 0;  % => mode = 0 
       xdot = []; % to prevent new warning msgs!!
   else
       % we are leaving "stuck" mode
       reset = [];  % no need to reset states...
       phi = sign(T_e);
       xdot = []; % to prevent new warning msgs!!
   end;
else
   error('bad value of mode')
end
