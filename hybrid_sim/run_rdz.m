% script to run rdz.m | rdz_seh.m to demonstrate state-event handling (SEH)
%
% JH Taylor - University of New Brunswick - 7 Mar 1996
%
% "rdz(x)" defined by sign(x) if |x| >= 1, else 0.
clear; close all;
t0 = 0.0; tf = 20.0; x0 = [ 2.5 ; 0.0 ]; dt = 0.05; 
timer = clock; 
[t1,x1] = trap_101('rdz_seh',t0,tf,x0,dt,eps,0); 
T_trap_101 = etime(clock,timer); % integration time with trapezoidal SEH
disp(['   TRAP_101 took ',num2str(T_trap_101),' seconds for two cycles']);
timer = clock;
[t2,x2] = ode45_101('rdz_seh',t0,tf,x0);%,29,1.e-6,0); % simulate 120 seconds!!
T_ode45_101 = etime(clock,timer); % integration time with state-event handling
disp(['   ODE45_101 took ',num2str(T_ode45_101),' seconds for two cycles']);
disp('Here is an example of the resulting waveforms . . .')
plot(t1,x1,t2,x2,'--') 
title('rdz\_seh: xdot(1) = x(2); xdot(2) = - rdz(x(1)) with Modes')
xlabel('time')
ylabel('x(t)')
text(2,2.1,'\_\_\_ = TRAP\_101')
text(7,-2.1,'\_ \_ = ODE45\_101')
disp('Note: the variable-step method ode45_101 is less smooth, due to taking')
disp('      large steps, but it may be substantially faster in some cases.')
disp('Strike any key to continue.')
pause   
disp('Here are the corresponding phase-plane portraits . . .');
figure(2); plot(x1(:,1),x1(:,2),x2(:,1),x2(:,2),'--');
title('rdz\_seh:  xdot(1) = x(2); xdot(2) = - rdz(x(1)) with Modes')
xlabel('x(1)')
ylabel('x(2)')
text(-1.10,0.15,'\_\_\_ = Algorithm TRAP\_101')
text(-1.13,-.15,'\_ \_ = Algorithm ODE45\_101') % IF tf = 20.0;
disp('Strike any key to continue.')
pause 
disp('The next step shows how much longer ODE45 takes - please be patient!!')
timer = clock;
[t3,x3] = ode45m('rdz',t0,tf,x0);
T_ode45 = etime(clock,timer); % integration time without seh
disp(['   so ODE45m took ',num2str(T_ode45),' seconds for two cycles!']);
disp('Note: ODE45m (m-file) is not compiled, to make the comparison fair . . .')
figure(3); plot(x1(:,1),x1(:,2),x2(:,1),x2(:,2),'--',x3(:,1),x3(:,2),'-.');
x_sw = [ 1 ; -sqrt(3) ];
hold on
plot(x_sw(1),x_sw(2),'o')
text(-1.10,0.25,['TRAP\_101 took ',num2str(T_trap_101),' seconds for two cycles']);
text(-1.10,0,['ODE45\_101 took ',num2str(T_ode45_101),' seconds for two cycles']);
text(-1.13,-.25,['ODE45m took ',num2str(T_ode45),' seconds for two cycles!'])
disp('At this scale results appear comparable . . .')
disp('Strike any key to continue.')
pause 
figure(4); plot(x1(:,1),x1(:,2),x2(:,1),x2(:,2),'--',x3(:,1),x3(:,2),'-.');
x_sw = [ 1 ; -sqrt(3) ];
hold on
plot(x_sw(1),x_sw(2),'o')
axis([1-1.e-03 1+3.5e-04 -1.7321 -1.7320]) % shows both switching points
text(1-7.e-04,-1.732047,'\_\_\_ = TRAP\_101')
text(1-5.2e-05,-1.732056,'o = Exact switching point')
text(1-5.2e-05,-1.732063,'      Pretty accurate, eh?')
text(1-4.e-04,-1.73208,'\_ . \_ . = ODE45M (m-file), first pass')
text(1-7.e-04,-1.732016,'ODE45M, another pass')
disp(' . . . but now you see the difference!')
disp(' ')
disp('End of demonstration.')
%% axis([1-5.e-07 1+5.e-07 -1.732051 -1.7320505]) % shows only SEH point
