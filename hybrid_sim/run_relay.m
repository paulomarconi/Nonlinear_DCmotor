% script to run relay_seh.m 
%
% JH Taylor - University of New Brunswick - 17 Feb 1996
%
clear; close all;
disp(' ')
disp('Here is the mode-based model for a periodically switching relay:')
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
type relay_seh
disp(' - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -')
disp(' ')
disp('Note: initializing the switching function correctly may be ')
disp('      complicated if you want to handle all cases carefully.')
disp(' ')
t0 = 0.0; tf = 20.0; x0 = [ .25 ; 0.0 ]; dt = 0.05; 
[t1,x1] = trap_101('relay_seh',t0,tf,x0,dt); 
disp('Here is an example of the resulting waveforms . . .')
plot(t1,x1)
title('relay\_seh:  xdot(1) = x(2); xdot(2) = - sign(x(1));')
xlabel('time')
ylabel('x(t)')
disp('Strike any key to continue.')
pause
disp('Here is the corresponding phase-plane portrait . . .');
figure(2); plot(x1(:,1),x1(:,2));
axis([-.3 .3 -.8 .8]);
title('relay\_seh:  xdot(1) = x(2); xdot(2) = - sign(x(1));')
xlabel('x(1)')
ylabel('x(2)')
text(-0.1,0.1,'\_\_\_ = Algorithm TRAP\_101')
if tf == 20.0, text(-0.1,-0.1,'Seven cycles of oscillation!'); end; 
disp('Strike any key to continue.')
pause
disp('Now we compare with the ode45m algorithm (ode45 m-file version) . . .');
hold on
%
[t2,x2] = ode45m('relay',t0,tf,x0);
plot(x2(:,1),x2(:,2),'--');
text(-.125,0.0,'\_ \_ = ODE45M (ode45 m-file version)')
disp('Note: the "jumps" are due to automatic use of a large step . . .');
disp('      we could use large steps with TRAP_101 too, it is exact for ');
disp('      these simple dynamics, but we wanted to show the waveform.');
%% Unfortunately RK45 in matlab 5 no longer allows this type of usage!
%% disp('Strike any key to continue.')
%% pause
%% disp('Now we compare with the SimuLink rk45 algorithm . . .');
%% hold on
%% %
%% [t3,x3] = rk45('relay_slink',tf,x0); removed/doesn't work...?
%% plot(x3(:,1),x3(:,2),'-.');
%% text(-.092,-.15,'\_ . \_ = SimuLink RK45')
disp('Strike any key to zoom into one switching point:')
pause
x_sw = [ 0 ; -sqrt(.5) ]; 
hold on
fudge = 5.e-4; f20 = fudge/20;
axis([-fudge fudge -sqrt(.5)-fudge -sqrt(.5)+fudge]);
plot(x_sw(1),x_sw(2),'o')
text(f20,-sqrt(.5)-f20,'o = Exact switching point')
text(f20,-sqrt(.5)-4*f20,'\_\_\_ = Algorithm TRAP\_101')
text(f20,-sqrt(.5)-7*f20,'\_ \_ = ODE45M (ode45 m-file version)')
disp('Strike any key to zoom much closer into the switching point:')
pause
figure(3); plot(x1(:,1),x1(:,2));
title('relay\_seh:  xdot(1) = x(2); xdot(2) = - sign(x(1));')
xlabel('x(1)')
ylabel('x(2)')
hold on
fudge = 5.e-11; f20 = fudge/20;
axis([-fudge fudge -sqrt(.5)-fudge -sqrt(.5)+fudge]);
plot(x_sw(1),x_sw(2),'o')
text(f20,-sqrt(.5)-f20,'o = Exact switching point,')
text(f20,-sqrt(.5)-4*f20,'      hit exactly 7 times . . . ')
text(f20,-sqrt(.5)-7*f20,'      Pretty accurate, eh?')
disp(['Note: this window is the exact point +/- ',num2str(fudge)])
