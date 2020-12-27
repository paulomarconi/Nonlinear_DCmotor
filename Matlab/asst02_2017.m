function xdot = asst02_2017(t,x)
global E_0 Tau_L0 T_Amb B_2C 

% motor parameters, Nachtigal, Table 16.5 p. 663

J_1 = 0.0035;   % in*oz*s^2/rad
B_1 = 0.064;    % in*oz*s/rad

% electrical/mechanical relations
K_E = 0.1785;       % back emf coefficient, e_m = K_E*omega_m (K_E=alpha*omega)
K_T = 141.6*K_E;    % torque coeffic., in English units K_T is not = K_E! (K_T=alpha*iA) 
R_A = 8.4;          % Ohms
L_A = 0.0084;       % H

% gear-train and load parameters
J_2 = 0.035;    % in*oz*s^2/rad % 10x motor J
B_2 = 2.64;     % in*oz*s/rad (viscous)
N = 8;          % motor/load gear ratio; omega_1 = N omega_2

% Thermal model parameters
R_TM = 2.2;     % Kelvin/Watt
C_TM = 9/R_TM;  % Watt-sec/Kelvin (-> 9 sec time constant - fast!)

Jeq = J_2+N*2*J_1;
Beq = B_2+N^2*B_1;
a = R_A/L_A;
b = K_E*N/L_A;
c = N*K_T/Jeq;
d = Beq/Jeq;
e = B_2C/Jeq;
f = R_A/C_TM;
g = 1/(C_TM*R_TM);

if t < 0.05
    e_i = 0;
else  
%     e_i = E_0;  
    e_i = E_0*sin(5*(2*pi)*(t - 0.05)); 
end
if t < 0.2
    Tau_L = 0;
else
    Tau_L = Tau_L0;
end

xdot(1) = -a*x(1)-b*x(2)+e_i/L_A;
xdot(2) = c*x(1)-d*x(2)-e*sign(x(2))-Tau_L/Jeq;
xdot(3) = f*x(1)^2-g*x(3)+g*T_Amb;
xdot = xdot(:); % force column vector
