clear; close all; clc;
global E_0 Tau_0 T_Amb B_2C ;

E_0 = 90; Tau_0 = 25; T_Amb = 18; B_2C = 140;
t0 = 0; tfinal = 1; stp = 0.0001;
x0 = [ 0; 0; 0 ]; %initial conditions

timer = clock;
[t1,x1] = ode45m('asst02_2017',t0,tfinal,x0,stp);
%[t1,x1] = eufix1('asst02_2017',[t0 tfinal],x0,stp);
Tsim1 = etime(clock,timer), % integration time 
Len1 = length(t1), % number of time-steps

figure

subplot(3,1,1)       
plot(t1,x1(:,1));
title('Motor+Load+thermal model simulation, B_{2C}=80')
xlabel('Time, t (sec)')
ylabel('i_A (amperes)')
text(0.2,1,['E_0 =',num2str(E_0,3)]);
%grid;

subplot(3,1,2)       
plot(t1,x1(:,2));
title('Motor+Load+thermal model simulation, B_{2C}=80')
xlabel('Time, t (sec)')
ylabel('\omega_2 , Motor angular velocity (rad/sec)')
text(0.2,24,['B_{2C} =',num2str(B_2C,3), '\tau_0 =',num2str(Tau_0,3)]);

subplot(3,1,3)      
plot(t1,x1(:,3));
title('Motor+Load+thermal model simulation, B_{2C}=80')
xlabel('Time, t (sec)')
ylabel('\theta_M (deg)')
text(0.2,5,['\Theta_A=',num2str(T_Amb,3) ' deg']);

