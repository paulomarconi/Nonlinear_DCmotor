J_1 = 0.0035; % in*oz*s^2/rad
B_1 = 0.064; % in*oz*s/rad
% electrical/mechanical relations
K_E = 0.1785; % back emf coefficient, e_m = K_E*omega_m
K_T = 141.6*K_E; % torque coeffic.; in English units K_T is not = K_E!
R_A = 8.4; % Ohms
L_A = 0.0084; % H
% gear-train and load parameters
J_2 = 0.035; % in*oz*s^2/rad % 10x motor J
B_2 = 2.64; % in*oz*s/rad (viscous)
N = 8; % motor/load gear ratio; omega_1 = N omega_2
% Thermal model parameters
R_TM = 2.2; % Kelvin/Watt
C_TM = 9/R_TM; % Watt-sec/Kelvin (-> 9 sec time constant - fast!)

E_0=120; Tau_0=80; T_Amb=25; B_2C=80;

Jeq=J_2+N^2*J_1;
Beq=B_2+N^2*B_1;

a=R_A/L_A;
b=K_T*N/L_A;
c=N*K_T;
d=Beq/Jeq;
e=B_2C/Jeq;
f=R_A/C_TM;
g=1/(C_TM*R_TM);



f = @(t,x) [-a*x(1)-b*x(2)+E_0/L_A;
            -c*x(1)-d*x(2)-e*sign(x(2))-Tau_0/Jeq;
            f*x(1)^2-g*x(3)+g*T_Amb;
            ];
[t,xa] = ode45(f,[0 10],[0 0 0]);

plot(t,xa(:,1))
title('x(t)')
xlabel('t'), ylabel('x')
pause;

plot(t,xa(:,2))
title('y(t)')
xlabel('t'), ylabel('y')
pause;

plot(t,xa(:,3))
title('z(t)')
xlabel('t'), ylabel('z')
