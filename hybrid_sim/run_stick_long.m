% script to run two models of a motor with stiction, one with
% state-event handling (SEH) and one without (ODE).
%
% JH Taylor - University of New Brunswick - 8 June 1996
%
clear; close all;
disp('  ')
disp('In this demonstration we use a second-order model that realistically')
disp('characterizes an electro-mechanical system with stiction and saturation. ')
disp('  ')
global v_in omega
v_in = 0.6; omega = 2.0;
t0 = 0.0; tf = 7.28; x0 = [ 0.0 ; 0.00 ]; dt = 0.05; 
timer = clock;
[t1,x1] = trap_101('stick_seh',t0,tf,x0,dt); 
T_trap_101 = etime(clock,timer); % integration time with trapezoidal SEH
disp([' TRAP_101 took ',num2str(T_trap_101),' seconds ; dt = ',num2str(dt)]);
timer = clock;
tol = 1.e-4; %% THIS IS 10x THE DEFAULT
[t2,x2] = ode45m('stick',t0,tf,x0,tol); 
T_ode45 = etime(clock,timer); % integration time with state-event handling
disp([' ODE45M took ',num2str(T_ode45),' seconds ; tol = ',num2str(tol)]);
plot(t1,x1(:,1),'--b',t1,x1(:,2),'-g',t2,x2(:,1),':r',t2,x2(:,2),'-.k') 
title('stiction model:  M*xdotdot = F_a - f_v*xdot - f_c*mode')
xlabel('time')
ylabel('x(t) = posn, xdot(t) = velo')
text(0.5,-10,['v_{in} = ',num2str(v_in)])
text(0.65,-13,['\omega = ',num2str(omega)])

hold on;
plot([1.00 1.7],[16.5 16.5],'--b',[1.9 2.6],[16.5 16.5],'-g');
text(1.75,16.5,',')
text(1.05,17.5,'(posn,'); text(2.05,17.5,'velo)'); 
text(2.8,16.5,[' via trap\_101, dt = ',num2str(dt)])

plot([3.30 4.0],[-17.5 -17.5],':r',[4.2 4.85],[-17.5 -17.5],'-.k');
text(4.1,-17.5,',')
text(3.35,-16.5,'(posn,'); text(4.25,-16.5,'velo)'); 
text(5.0,-17.5,[' via ode45m, tol = ',num2str(tol)])
%% axis([0 6.4 -18 18])
disp('Note: this comparison uses ODE45M (the matlab m-file), which is uncompiled,')
disp('      to obtain a meaningful CPU time comparison.  Also, note that to keep')
disp(['      the run-time manageable we set tol = ',num2str(tol),' which gave rise to the'])
disp('      rough chattering when the velocity should be == 0.')
