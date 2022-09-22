clear all; 
close all; 
clc

format long

n = 60*60;  % Number of second time points 

% Heater steps
Q1 = ones(n,1)*30;
Q2 = ones(n,1)*30;

Q_matrix=[Q1 Q2];

% Initial temperature
T0 = 23.0 + 273.15;

% Store temperature results
T1 = ones(n,1)*T0;
T2 = ones(n,1)*T0;

time = linspace(0,n,n); % Time vector

for i = 2:n
    % initial condition for next step
    x0 = [T1(i-1),T2(i-1)];
    % time interval for next step
    tm = [time(i-1),time(i)];
    % Integrate ODE for 1 sec each loop
    z = ode45(@(t,x)heat2(t,x,Q1(i-1),Q2(i-1)),tm,x0); %Q1 and Q2 at i-1 cuz we are wanting things to be calculated for 1 second
    % record T1 and T2 at end of simulation
    T1(i) = z.y(1,end);
    T2(i) = z.y(2,end);
end

T_predicted_matrix=[T1 T2];

% Plot results
figure(1)

subplot(2,2,1)
plot(time/60.0,Q1,'b--','LineWidth',2)
ylabel('Heater Output')
legend('Q_1')
ylim([0,50])
xlabel('Time (min)')

subplot(2,2,2)
plot(time/60.0,T1-273.15,'b-','LineWidth',2)
hold on
plot(time(end)/60.0,T1(end)-273.15,'bo')
legend('T1 predicted')

subplot(2,2,3)
plot(time/60.0,Q2,'r--','LineWidth',2)
ylabel('Heater Output')
legend('Q_2')
ylim([0,50])
xlabel('Time (min)')

subplot(2,2,4)
plot(time/60.0,T2-273.15,'r-','LineWidth',2)
ylabel('Temperature (degC)');
legend('T2 predicted')
ylim([20 40]);

sgtitle('Temperature simulation for Q1=30 and Q2=30 for t=60s')

fprintf('Error between last steady states value in T1 %.16f\n', T1(end)-T1(end-1));
fprintf('Error between last steady states value in T2 %0.16f\n', T2(end)-T2(end-1));

% save as heat2.m
% define energy balance model
function dTdt = heat2(t,x,Q1,Q2)
    % Parameters
    Ta = 23 + 273.15;   % K
    U = 10.0;           % W/m^2-K
    m = 4.0/1000.0;     % kg
    Cp = 0.5 * 1000.0;  % J/kg-K    
    A = 10.0 / 100.0^2; % Area in m^2
    As = 2.0 / 100.0^2; % Area in m^2
    alpha1 = 0.0100;    % W / % heater 1
    alpha2 = 0.0075;    % W / % heater 2
    eps = 0.9;          % Emissivity
    sigma = 5.67e-8;    % Stefan-Boltzman

    % Temperature States
    T1 = x(1);
    T2 = x(2);

    % Heat Transfer Exchange Between 1 and 2
    conv12 = U*As*(T2-T1);
    rad12  = eps*sigma*As * (T2^4 - T1^4);

    % Nonlinear Energy Balances
    dT1dt = (1.0/(m*Cp))*(U*A*(Ta-T1) ...
            + eps * sigma * A * (Ta^4 - T1^4) ...
            + conv12 + rad12 ...
            + alpha1*Q1);
    dT2dt = (1.0/(m*Cp))*(U*A*(Ta-T2) ...
            + eps * sigma * A * (Ta^4 - T2^4) ...
            - conv12 - rad12 ...
            + alpha2*Q2);
    dTdt = [dT1dt,dT2dt]';
end
