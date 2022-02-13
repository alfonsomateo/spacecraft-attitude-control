% Spacecraft Attitude Dynamics and Control

clear all;
clc;

% Simulation parameters

mu = 398600.441; % km^3/s^2
R_E = 6378.136; % km
h = 700; % km
n = sqrt(mu/(R_E+h)^3);

J = [124.531, 0, 0;
    0, 124.586, 0;
    0, 0, 1.704]; % kgm^2
Td = [1e-4; 1e-4; 1e-4]; % Nm
ts = 0.1; % s
tf = 1500; % s
initial_att = deg2rad([30, 30, 30]);


att_commands_t = [0,100,500.1,900.1,1500.1];
att_commands_Euler = deg2rad([0, 0, 0;
    60, 60, 60;
    -60, -60, -60;
    0, 0, 0]);
sz = size(att_commands_Euler);
att_commands_quaternions = zeros(sz(1),4);
for i = 1:sz(1)
    quat = eul2quat(att_commands_Euler(i,:),'ZYX');
    att_commands_quaternions(i,:) = [quat(2:4) quat(1)];
end

% Select between Euler angles and quaternions

att_repr = 'Euler'; % select between using Euler angles or quaternions

% Control parameters

% Kp = 0.02*J;
% Kd = 0.5*J;
Kp = 0.02*eye(3);
Kd = 0.5*eye(3);

control_type = 'PD';

% Integration parameters

dt = ts/10;

% Save all simulation parameters into a struct

sim = struct('n',n,'J',J,'Td',Td,'ts',ts,'tf',tf,...
    'initial_att',initial_att,'att_commands_t',att_commands_t,...
    'att_commands_Euler',att_commands_Euler,'att_commands_quaternions',att_commands_quaternions,...
    'att_repr',att_repr,'Kp',Kp,'Kd',Kd,'control_type',control_type,'dt',dt);

%% Run simulations

% PD Euler
sim1 = sim;
sim1.att_repr = 'Euler';
sim1.control_type = 'PD';
sim1.Kp = 0.02*J/2;
sim1.Kd = 0.5*J/2;

% PD quaternions
sim2 = sim;
sim2.att_repr = 'quaternions';
sim2.control_type = 'PD';
sim2.Kp = 0.02*J;
sim2.Kd = 0.5*J;

% NDI Euler
sim3 = sim;
sim3.att_repr = 'Euler';
sim3.control_type = 'NDI';
sim3.Kp = 0.02*eye(3)/2;
sim3.Kd = 0.5*eye(3)/2;

% NDI quaternions
sim4 = sim;
sim4.att_repr = 'quaternions';
sim4.control_type = 'NDI';
sim4.Kp = 0.02*eye(3);
sim4.Kd = 0.5*eye(3);

% NDI TSSP Euler
sim5 = sim;
sim5.att_repr = 'Euler';
sim5.control_type = 'NDI_tss';
sim5.Kp = 0.02*eye(3)*2;
sim5.Kd = 0.5*eye(3)/2;

% NDI TSSP quaternions
sim6 = sim;
sim6.att_repr = 'quaternions';
sim6.control_type = 'NDI_tss';
sim6.Kp = 0.02*eye(3)*2;
sim6.Kd = 0.5*eye(3)/2;

% INDI TSSP Euler
sim7 = sim;
sim7.att_repr = 'Euler';
sim7.control_type = 'INDI';
sim7.Kp = 0.02*eye(3)*2;
sim7.Kd = 0.5*eye(3)/2;

% INDI TSSP quaternions
sim8 = sim;
sim8.att_repr = 'quaternions';
sim8.control_type = 'INDI';
sim8.Kp = 0.02*eye(3)*2;
sim8.Kd = 0.5*eye(3)/2;

sims = [sim1,sim2,sim3,sim4,sim5,sim6,sim7,sim8];

parfor i = 1:8
    S(i) = run_dynamics(sims(i));
end

%% Plots

% Display simulation attitude results
att_ref_Euler = zeros(15001,3);
att_ref_quaternions = zeros(15001,4);
time_ref = 0:ts:tf;
for i = 1:15001
    [att_c,att_dot_c] = comanded_attitude(time_ref(i),sim1);
    att_ref_Euler(i,:) = att_c';
    [att_c,att_dot_c] = comanded_attitude(time_ref(i),sim2);
    att_ref_quaternions(i,:) = att_c';
end
figure(1001)
subplot(3,1,1)
plot(time_ref,rad2deg(att_ref_Euler(:,1)),'Color',[0 0.4470 0.7410],'LineWidth',1.2)
ylim([-90,90])
xlabel('Time [s]');ylabel('\theta_{1c} [deg]');grid
set(gca,'FontSize',12)
subplot(3,1,2)
plot(time_ref,rad2deg(att_ref_Euler(:,2)),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.2)
ylim([-90,90])
xlabel('Time [s]');ylabel('\theta_{2c} [deg]');grid
set(gca,'FontSize',12)
subplot(3,1,3)
plot(time_ref,rad2deg(att_ref_Euler(:,3)),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2)
ylim([-90,90])
xlabel('Time [s]');ylabel('\theta_{3c} [deg]');grid
set(gca,'FontSize',12)
sgtitle('Reference commands')

figure(1002)
subplot(4,1,1)
plot(time_ref,att_ref_quaternions(:,1),'Color',[0 0.4470 0.7410],'LineWidth',1.2)
ylim([-1,1])
xlabel('Time [s]');ylabel('q_{1c} [-]');grid
set(gca,'FontSize',12)
subplot(4,1,2)
plot(time_ref,att_ref_quaternions(:,2),'Color',[0.8500 0.3250 0.0980],'LineWidth',1.2)
ylim([-1,1])
xlabel('Time [s]');ylabel('q_{1c} [-]');grid
set(gca,'FontSize',12)
subplot(4,1,3)
plot(time_ref,att_ref_quaternions(:,3),'Color',[0.9290 0.6940 0.1250],'LineWidth',1.2)
ylim([-1,1])
xlabel('Time [s]');ylabel('q_{1c} [-]');grid
set(gca,'FontSize',12)
subplot(4,1,4)
plot(time_ref,att_ref_quaternions(:,4),'Color',[0.4940 0.1840 0.5560],'LineWidth',1.2)
ylim([-1,1])
xlabel('Time [s]');ylabel('q_{4c} [-]');grid
set(gca,'FontSize',12)
sgtitle('Reference commands')

titles = ["PD Euler",'PD quaternions','NDI Euler','NDI quaternions','NDI TSSP Euler','NDI TSSP quaternions','INDI TSSP Euler','INDI TSSP quaternions'];
for i = 1:8
    figure(i)
    subplot(1,2,1)
    plot(S(i).t,rad2deg(S(i).Euler(:,1)),'-','LineWidth',1.2)
    hold on
    plot(S(i).t,rad2deg(S(i).Euler(:,2)),'--','LineWidth',1.2)
    plot(S(i).t,rad2deg(S(i).Euler(:,3)),'-.','LineWidth',1.2)
    h = plot(time_ref,rad2deg(att_ref_Euler(:,1)),'k');
    hold off
    ylim([-90,90])
    legend('\theta_1','\theta_2','\theta_3','\theta_{ic}')
    xlabel('Time [s]');ylabel('Attitude Euler angles [deg]');grid
    uistack(h,'bottom')
    set(gca,'FontSize',12)
    subplot(1,2,2)
    plot(S(i).t,S(i).quaternions(:,1),'-','LineWidth',1.2)
    hold on
    plot(S(i).t,S(i).quaternions(:,2),'--','LineWidth',1.2)
    plot(S(i).t,S(i).quaternions(:,3),'-.','LineWidth',1.2)
    plot(S(i).t,S(i).quaternions(:,4),':','LineWidth',1.2)
    h = plot(time_ref,att_ref_quaternions,'k');
    hold off
    ylim([-1,1])
    legend('q_1','q_2','q_3','q_4','q_{ic}')
    xlabel('Time [s]');ylabel('Attitude quaternions [-]');grid
    uistack(h,'bottom')
    set(gca,'FontSize',12)
    sgtitle(titles(i))
end

figure(9)
plot(S(1).t,rad2deg(S(1).theta),'LineWidth',1.2)
hold on
for i = 2:8
    plot(S(i).t,rad2deg(S(i).theta),'LineWidth',1.2)
end
hold off
xlim([500.2,600]);ylim([0,180])
legend(titles)
xlabel('Time [s]');ylabel('\theta [deg]');grid
set(gca,'FontSize',12)

figure(10)
plot(S(1).t,S(1).Tc_norm,'LineWidth',1.2)
hold on
for i = 2:8
    plot(S(i).t,S(i).Tc_norm,'LineWidth',1.2)
end
hold off
xlim([500.2,550])
legend(titles)
xlabel('Time [s]');ylabel('T_c [Nm]');grid
set(gca,'FontSize',12)

titles2 = ["PD Euler",'PD quaternions','NDI Euler','NDI quaternions','NDI & INDI TSSP Euler','NDI & INDI TSSP quaternions'];

figure(11)
plot(S(1).t,rad2deg(S(1).theta),'LineWidth',1.2)
hold on
for i = 2:6
    plot(S(i).t,rad2deg(S(i).theta),'LineWidth',1.2)
end
hold off
xlim([500.2,600]);ylim([0,180])
legend(titles2)
xlabel('Time [s]');ylabel('\theta [deg]');grid
set(gca,'FontSize',12)

figure(12)
plot(S(1).t,S(1).Tc_norm,'LineWidth',1.2)
hold on
for i = 2:6
    plot(S(i).t,S(i).Tc_norm,'LineWidth',1.2)
end
hold off
xlim([500.2,550])
legend(titles2)
xlabel('Time [s]');ylabel('T_c [Nm]');grid
set(gca,'FontSize',12)

figure(13)
plot(S(1).t,S(1).omega_norm,'LineWidth',1.2)
hold on
for i = 2:6
    plot(S(i).t,S(i).omega_norm,'LineWidth',1.2)
end
hold off
xlim([500.2,600])
legend(titles2)
xlabel('Time [s]');ylabel('\omega [rad/s]');grid
set(gca,'FontSize',12)

