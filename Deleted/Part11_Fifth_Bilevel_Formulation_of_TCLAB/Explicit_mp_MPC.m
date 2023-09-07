clc;
clear;

%% loaded all the necessary files

Matrices_defined_for_Twolevel;
load('Q1_30_Q2_30.mat') %steady state data
T1_measured_initial=T1(end,1); %steady state value
T2_measured_initial=T2(end,1); %steady state value
ynominal=[314.4435277118964	311.146230339323];

tic

%% Defined all initial conditions needed

%Initially we kept heaters on for 30% and allowed the system to reach steady state
u0 = [0 0];%input  
u = u0;

%Initial values of states
%Taking a row vector of xSSp
xSSp = ss1.Report.Parameters.X0';%pseudo states from state space model
%Initial temperature values
y0 = [ss1.C*xSSp']' + ynominal;
% y0(1)=314.4434;
% y0(2)=311.1456;
% y0(1)=310;
% y0(2)=315;
% y0(1)= 320;
% y0(2)= 315;
y0(1)=300;
y0(2)=300;
%Giving the steady state temperatures as set point 
ysp = [T1_measured_initial T2_measured_initial];

%% Defined all variable matrices that will change after each sample time

%Running simulation for 150 seconds
timeTotal = 1000;% total time(sec)

%Defining initial matrices and conditions 
t0 = 0;
ymat = y0; % Matrix for storing values of controlled Output after each sample time %
umat = u0;% Matrix for storing control input after each sample time
xSSpmat = xSSp; %Matrix for storing states after each sample time 
tmat = t0;
yspmat = ysp;%Matrix for storing values of set points after each sample time

%% Running for loop , taking sample time as one second

t_int = 1; % sec

% %   Discrete-time identified state-space model:
% %     x(t+Ts) = A x(t) + B u(t) + K e(t)
% %        y(t) = C x(t) + D u(t) + e(t)

% % problem.namesThita
tt = [];
for i=0:t_int:timeTotal-t_int
    tspan=[i i+t_int];

    %Defined parameter using the definition in problem.namesThita
    theta1 = [xSSp y0-ynominal ysp-ynominal]; % 
   
   [nCR1,uaux1] = PointLocation(Solution_outer,theta1');
   %[uaux,fval,exitflag] = cplexqp(2*problem.Q, problem.Ht*theta'+problem.c, problem.A, problem.b+problem.F*theta');
    %uaux1 = 0;
    theta2=[xSSp uaux1 y0-ynominal ysp-ynominal];
    [ncR2,uaux2]= PointLocation(Solution_inner,theta2');
    uaux=[uaux1;uaux2];
    
    uaux = real(uaux); %to discard 0.00001i due to numerical errors
    if(isempty(uaux))
        disp('uaux empty break')
        break
    end
    
    u = uaux(1:2); % inputs at the current time step
    umat = [umat;u'];
    
    xSSp = ss1.A*xSSpmat(end,:)' + ss1.B(:,2)*uaux2 + ss1.B(:,1)*uaux1;%*(umat(end,:))'); %no disturbance
    ySSp = (ss1.C*xSSpmat(end,:)');
    
    %calculating output for optimal control inputs obtained
    [t,y] = ode45(@(t,x)heat2Model(t,x,u'+unominal),tspan,y0);
    y0 = y(end,:);
    
    xSSpmat = [xSSpmat;xSSp'];
    ymat = [ymat; y0];
    yspmat = [yspmat;ysp];
    
    tmat = [tmat; i+t_int];
    tt = [tt;nCR1 ncR2];
    xSSp=xSSp';
end

umat = umat + unominal;
simElapsedTime = toc;

%% To check whether MPC tuning is successful 

if ((ymat(end,1)-yspmat(end,1))==0 && (ymat(end,2)-yspmat(end,2))==0)
    disp("Tuning Successful")
else
    disp("tuning not successful")
end
%% plots

figure(2)

%Plots for temperatures
subplot(2,2,1)
plot(tmat,ymat(:,1),'k','linewidth',3)
hold all
plot(tmat,yspmat(:,1),'k','linestyle','--','linewidth',3)
hold off
% axis([0 timeTotal 305 330])
xlabel('time(s)','FontSize',20)
ylabel('T1','FontSize',20)
%ylim([300,330])
set(gca,'Fontsize',20)

subplot(2,2,2)
plot(tmat,ymat(:,2),'k','linewidth',3)
hold all
plot(tmat,yspmat(:,2),'k','linestyle','--','linewidth',3)
hold off
% axis([0 timeTotal 305 330])
xlabel('time(s)','FontSize',20)
ylabel('T2','FontSize',20)
%ylim([300,330])
set(gca,'Fontsize',20)

subplot(2,2,3)
plot(tmat,umat(:,1),'k','linewidth',3)
% axis([0 timeTotal 0 umax(1)])
xlabel('time(s)','FontSize',20)
ylabel('u1','FontSize',20)
set(gca,'Fontsize',20)

subplot(2,2,4)
plot(tmat,umat(:,2),'k','linewidth',3)
% axis([0 timeTotal 0 umax(2)])
xlabel('time(s)','FontSize',20)
ylabel('u2','FontSize',20)
set(gca,'Fontsize',20)

%sgtitle('Multiparametric MPC: Tuning Parameters are QR=diag[9300,9300],  R1=diag[0.5,0.5],  OH=2,  NC=2')
% sgtitle('MPC: Tuning Parameters are QR=diag[5000,5000],  R1=diag[0.5,0.5],  OH=2,  NC=2')