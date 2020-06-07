% script to simulate Gibson's roll control problem using ode45_101
% JH Taylor - 24 March 1999

%% uses model with SI units
%% NB: figure 2 is not published; it shows the 
%% exact capture of the switching times
clear; close all;
% initial conditions etc.
x0 = [ -0.1621 -5.0781 0.0815 1.184]; t0 = 0.0; tf = 2.1; % ~6 cycles
[t,x] = ode45_101('gibson_si',t0,tf,x0);
figure(1)
subplot(211),plot(t,x(:,1),'-'),axis([0 tf -.4 .4]); 
xlabel('time'); ylabel('\phi');
title('Roll control problem from Gibson -- \phi(t)')
hold on                                             
plot([ 1.4 tf ],[ 0.133 0.133 ],'--')                
plot([ 0 tf ],[ 0 0 ],':')                
plot([ 1.4 tf ],[ -0.133 -0.133 ],'--')
text(1.6,0.2,'\phi = 0.133')
plot([1.57 1.52],[ .2 .135])
plot([ 1.6614 1.6614 ],[ -0.2 0.03 ],'-.') % taken from actual switching times
plot([ 1.9335 1.9335 ],[ -0.2 0.03 ],'-.')
text(1.54,-.24,'t = 1.6614')
text(1.8,-.24,'t = 1.9335')
T = 1.9335 - 1.6614; omega = 2*pi/T; text(1.64,-0.32,['\omega = ',num2str(omega)])
subplot(212),plot(x(:,1),x(:,2),'-');
xlabel('\phi'); ylabel('\phi - dot ');
title('Roll control problem from Gibson -- \phi-dot vs \phi')
%
figure(2)
plot(t,2000*x(:,3),'-'); %% 2000*x(3) = relay input
hold on;
h = 22.24; %% thruster "relay" constant (Newtons)
plot([ 0 tf ],[ h h ],'--')
plot([ 0 tf ],[ -h -h ],'--')
plot([ 0.19137 0.19137],[-50 50],'--') %% these are the exact switching times
plot(0.19137,-h,'o');text(0.2,-20,'First switching point');
plot([ 0.38069 0.38069],[-50 50],'--') %% reported by ode45_101
plot(0.38069, h,'o');text(0.39,20,'Second switching point');
plot([ 0.54358 0.54358],[-50 50],'--')
plot(0.54358,-h,'o');text(0.56,-20,'Third switching point');
plot([ 0.69377 0.69377],[-50 50],'--')
plot([ 0.84608 0.84608],[-50 50],'--')
title('Roll control problem from Gibson -- 2000*x(3) = relay input')
xlabel('time'); ylabel('e(t) = 2000*x_3(t)');
axis([ 0.15 0.5 -1.2*h 1.2*h ])
%
% figure(1);
% print -deps 2409fg08
