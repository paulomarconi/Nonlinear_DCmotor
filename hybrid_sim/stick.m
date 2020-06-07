function xdot = stiction(t,x)
% model of an electro-mechanical system with stiction and saturation
%
% simple "relay" version of Taylor & Strobel, 1985 ACC paper.
%                                           JH Taylor, 16 June 1996.
%
global v_in delta omega
%
m_1 = 5.0;     % Nm/v
delta = 0.5;   % v
m_2 = 1.0;     % Nm/v
f_v = 0.1;     % Nms/rad
f_c = 1.0;     % Nm
MoI = 0.01;    % kg-m^2
% first, prelim calcs: 
% saturate the input voltage -> electrical torque:
volts = v_in*sin(omega*t);
av = abs(volts); sv = sign(volts);
if av < delta, T_e = m_1*volts;
else T_e = (m_1*delta + m_2*(av-delta))*sv; end 
%
xdot(1) = x(2);
xdot(2) = ( T_e - f_v*x(2) - f_c*sign(x(2)) )/MoI;
