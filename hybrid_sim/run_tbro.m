% script to demonstrate ode45_101 (ode45 integration with state-event
% handling including state RESET and SIMULTANEOUS state events)
%
% JH Taylor - 4 March 1996 - revised 2 February 2005
%
%% NB: the switching intervals become shorter and shorter, approaching
%% zero as t -> 9*sqrt(2*x0); e.g. for x0(1) = 0.25 MAX tf = 6.363961031!
%% (based on a reset of 0.8 each switch, see twin_ball_reset.m)

clear; close all;
disp('  ');
disp('Note: this demonstration does not include a standard ode45 simulation,');
disp('      since it cannot handle state variable reset. ');
disp('  ');
disp('This will take a few seconds, there are a lot of switchings ... ');
disp('  ');
t0 = 0.0; tf = 6.3639;
x0 = [ .25; 0.0; -.25; 0.0 ];
% calling seq: ode45_101(dyfun, t0, tf, y0, tend, tol, trace)
[t,x] = ode45_101('twin_ball_reset',t0,tf,x0); 
subplot(211),plot(t,x(:,1:2),t,x(:,3:4),'--');
title('twin\_ball\_reset (two "relays" with state reset); simultaneous SEs');
xlabel('time');
ylabel('x(t)');
subplot(212),plot(x(:,1),x(:,2),x(:,3),x(:,4),'--');
title('ODE45\_101 = 3-level mode state-event handler');
xlabel('x(1), x(3)');
ylabel('x(2), x(4)');
hold on
subplot(211),text(2.2,-.75,['date = ',date]);
subplot(212),plot(.25,0,'o');
  text(0.21,0.20,['x0(1) = ',num2str(x0(1),6)]);
subplot(212),plot(-.25,0,'o');
  text(-.29,-.20,['x0(3) = ',num2str(x0(3),6)]);
subplot(212),axis([-.3,.3,-1,1]);
disp('  ');
disp('Note: the switching rate becomes infinite as t -> 6.363961!');
disp('   There are 52 switchings in the interval [0 6.3639] . . .');
disp('  ');
disp('Strike any key to continue.');
pause
disp('  ');
disp('   Here is a more detailed view of the end-game switching:');
figure(2);
subplot(211),plot(t,[500*x(:,1) x(:,2)],t,[500*x(:,3) x(:,4)],'--');
subplot(211),title('twin\_ball\_reset, last 30 simultaneous SEs');
subplot(211),axis([6.32 6.364 -0.007 0.007])
subplot(212),plot(x(:,1),x(:,2),x(:,3),x(:,4),'--');
subplot(212),axis([-.001 .001 -0.05 0.05])
disp('  ');
disp('End of demonstration')
