clear variables; close all; clc;
global E_0 Tau_L0 T_Amb B_2C;

E_0 = 120; % [V]        120
Tau_L0 = 80; % [N.m]    80
T_Amb = 18; % [deg]     18
B_2C = 300; % [N]        80/300

t0 = 0; tfinal = 0.3; step = 1e-4;
x0 = [0; 0; 0]; % initial conditions

input_type = 0; % 0=constant, 1=sinusoidal
%% ode45 vs ode45m vs eufix1

timer = clock; 
[t1,x1] = ode45('asst02_2017',[t0, tfinal],x0);
% [t1,x1] = ode45m('asst02_2017',t0,tfinal,x0,step);
Tsim1 = etime(clock,timer);  % integration time 
Len1 = length(t1);           % number of time-steps 

timer = clock;
[t2,x2] = ode45m('asst02_2017',t0,tfinal,x0,step);
Tsim2 = etime(clock,timer);  % integration time 
Len2 = length(t2);           % number of time-steps

timer = clock;
[t3,x3] = eufix1('asst02_2017',[t0 tfinal],x0,step);
Tsim3 = etime(clock,timer);  % integration time 
Len3 = length(t3);           % number of time-steps

%% Relative error

% relative error at max current: ode45 vs eufix1
max_iA_ode45 = max(x1(:,1)); 
max_iA_eufix1 = max(x3(:,1)); 
max_iA_error = 100*abs( (max_iA_ode45-max_iA_eufix1)/max_iA_ode45 ) ; 

% relative error at max angular velocity: ode45 vs eufix1
max_omega2_ode45 = max(x1(:,2)); 
max_omega2_eufix1 = max(x3(:,2)); 
max_omega2_error = 100*abs( (max_omega2_ode45-max_omega2_eufix1)/max_omega2_ode45 ); 

%% Plotting
if input_type == 0
    %% Constant input e_i=E0
    figure;
        subplot(3,1,1);       
        plot(t1,x1(:,1),t2,x2(:,1),'--',t3,x3(:,1),'-.','LineWidth',1.5);
        title(['Nonlinear DC motor with thermal model, $B_{2C}=$',num2str(B_2C)],'Interpreter','Latex');
        ylabel('$i_A$ [A]','Interpreter','Latex');
        legend(['ode45: ',num2str(Tsim1),' [s]'],['ode45m: ',num2str(Tsim2),' [s]'],['eufix1: ',num2str(Tsim3),' [s]']);
        grid on;

        subplot(3,1,2);       
        plot(t1,x1(:,2),t2,x2(:,2),'--',t3,x3(:,2),':','LineWidth',1.5);
        ylabel('$\omega_2$ [rad/s]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

        subplot(3,1,3);      
        plot(t1,x1(:,3),t2,x2(:,3),'--',t3,x3(:,3),':','LineWidth',1.5);
        xlabel('Time [s]','Interpreter', 'Latex');
        ylabel('$\theta_M$ [deg]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

    % print('../asst02_2017/E0_ode45-ode45m-eufix1_1e-3.png','-dpng','-r300'); % Save as PNG with 300 DPI            

    figure;
        subplot(2,1,1);
        plot(t1,x1(:,1),t2,x2(:,1),'--',t3,x3(:,1),'-.','LineWidth',1.5);
        title(['Nonlinear DC motor with thermal model, $B_{2C}=$',num2str(B_2C)],'Interpreter','Latex');
        ylabel('$i_A$ [A]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1');
        axis([0.05 0.07 -inf inf]);
        text(0.058,5.5,['Relative error at max $i_{A}$=',num2str(max_iA_error),' $\%$'],'Interpreter','Latex');
        grid on;

        subplot(2,1,2);       
        plot(t1,x1(:,2),t2,x2(:,2),'--',t3,x3(:,2),':','LineWidth',1.5);
        ylabel('$\omega_2$ [rad/s]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        axis([0.05 0.07 -inf inf]);
        text(0.058,40,['Relative error at max $\omega_{2}$=',num2str(max_omega2_error),' $\%$'],'Interpreter','Latex');
        grid on;

    % print('../asst02_2017/E0_ode45-ode45m-eufix1_1e-3_zoom.png','-dpng','-r300'); % Save as PNG with 300 DPI     

    figure;
        plot(t1,x1(:,3),t2,x2(:,3),'--',t3,x3(:,3),':','LineWidth',1.5);
        title('Motor temperature $\theta_M$ over $80[s]$','Interpreter','Latex');
        xlabel('Time [s]','Interpreter', 'Latex');
        ylabel('$\theta_M$ [deg]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

    % print('../asst02_2017/thetaM_ode45-ode45m-eufix1_1e-3.png','-dpng','-r300'); % Save as PNG with 300 DPI     

elseif input_type == 1
    %% Sinusoidal input e_i

    figure;
        subplot(3,1,1);       
        plot(t1,x1(:,1),t2,x2(:,1),'--',t3,x3(:,1),'-.','LineWidth',1.5);
        title(['Nonlinear DC motor with thermal model, $B_{2C}=$',num2str(B_2C)],'Interpreter','Latex');
        ylabel('$i_A$ [A]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

        subplot(3,1,2);       
        plot(t1,x1(:,2),t2,x2(:,2),'--',t3,x3(:,2),':','LineWidth',1.5);
        ylabel('$\omega_2$ [rad/s]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

        subplot(3,1,3);      
        plot(t1,x1(:,3),t2,x2(:,3),'--',t3,x3(:,3),':','LineWidth',1.5);
        xlabel('Time [s]','Interpreter', 'Latex');
        ylabel('$\theta_M$ [deg]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

    % print('../asst02_2017/sinE0_ode45-ode45m-eufix1_1e-4.png', '-dpng', '-r300'); % Save as PNG with 300 DPI            

    figure;
        subplot(3,1,1);
        plot(t1,x1(:,2),t2,x2(:,2),'--',t3,x3(:,2),'-.','LineWidth',1.5);
        title(['Stiction behaviour on $\omega_2$, $B_{2C}=$',num2str(B_2C)],'Interpreter','Latex');
        ylabel('$\omega_2$ [rad/s]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        axis([0.148 0.157 -0.6 0.4]);
        grid on;

        subplot(3,1,2);       
        plot(t1,x1(:,2),t2,x2(:,2),'--',t3,x3(:,2),':','LineWidth',1.5);
        ylabel('$\omega_2$ [rad/s]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        axis([0.148 0.157 -0.015 0.010]);
        grid on;

        subplot(3,1,3);
        plot(t1,x1(:,2),t2,x2(:,2),'--',t3,x3(:,2),'-.','LineWidth',1.5);
        xlabel('Time [s]','Interpreter', 'Latex');
        ylabel('$\omega_2$ [rad/s]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        axis([0.148 0.157 -11e-5 5e-5]);
        grid on;

    % print('../asst02_2017/sinE0_ode45-ode45m-eufix1_1e-4_zoom.png', '-dpng', '-r300'); % Save as PNG with 300 DPI     

    figure;
        plot(t1,x1(:,3),t2,x2(:,3),'--',t3,x3(:,3),':','LineWidth',1.5);
        title('Motor temperature $\theta_M$ over $80[s]$','Interpreter','Latex');
        xlabel('Time [s]','Interpreter', 'Latex');
        ylabel('$\theta_M$ [deg]','Interpreter','Latex');
        legend('ode45','ode45m','eufix1','Location','southeast');
        grid on;

    % print('../asst02_2017/sinE0_thetaM_ode45-ode45m-eufix1_1e-4.png','-dpng','-r300'); % Save as PNG with 300 DPI     
end