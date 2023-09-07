clc;
clear;

%% loaded all the necessary files

Matrices_defined_for_comparing_MPC_and_mpMPC;
load('Q1_30_Q2_30.mat')
T1_measured_initial=T1(end,1);
T2_measured_initial=T2(end,1);

% %Reformulation of the problem
problem = ss2qp_yalmip(mpc_crystal,model_crystal);
problem.Q = problem.Q/2;

tic

%% Defined all initial conditions needed

%Initially we kept heaters on for 30% and allowed the system to reach
%steady state
u0 = 0*ones(1,2);%input  
u = u0;

%Initial values of states
%Taking a row vector
xSSp = ss1.Report.Parameters.X0'*0.92; %pseudo states from state space model
%Initial temperature values
y0 = [ss1.C*xSSp']';
% y0(1)=320;
% y0(2)=300;
%Giving the steady state temperatures as set point 
ysp = [T1_measured_initial T2_measured_initial];%*1.05; %SP change after 2 s

%% Defined all variable matrices that will change after each sample time

%Running simulation for 150 seconds
timeTotal = 200;% total time(sec)

%Defining initial matrices and conditions 
t0 = 0;
ymat = y0; % Matrix for storing values of controlled Output after each sample time %
umat = u0;% Matrix for storing control input after each sample time
xSSpmat = xSSp; %Matrix for storing states after each sample time 
tmat = t0;
yspmat = ysp;%Matrix for storing values of set points after each sample time

%% Running for loop , taking sample time as one second

t_int = 1; % sec

% from output
% % statespaceTanks =
% %   Discrete-time identified state-space model:
% %     x(t+Ts) = A x(t) + B u(t) + K e(t)
% %        y(t) = C x(t) + D u(t) + e(t)

% % problem.namesThita
% % 
% % ans =
% % 

for i=0:t_int:timeTotal-t_int
    tspan=[i i+t_int];%     if i == 20
%         %Giving step change to temperature values
%         ysp = 1.5*ysp;   
%     end
    %Defined parameter using the definition in problem.namesThita
    theta = [xSSp ymat(end,:) ysp umat(end,:)];
    
    [nCR,uaux,fval1] = PointLocation(Solution,theta');
    %[uaux,fval,exitflag] = cplexqp(2*problem.Q, problem.Ht*theta'+problem.c, problem.A, problem.b+problem.F*theta');

    uaux = real(uaux); %to discard 0.00001i due to numerical errors
    if(isempty(uaux))
        disp('uaux empty break')
        break
    end
    
    u = uaux(1:2); % inputs at the current time step
    umat = [umat;u'];
    
    xSSp = (ss1.A*xSSpmat(end,:)' + ss1.B*umat(end,:)')'; %no disturbance
    ySSp = ss1.C*xSSpmat(end,:)';
    
    xSSpmat = [xSSpmat;xSSp];
    ymat = [ymat; ySSp'];
    yspmat = [yspmat;ysp];
    
    tmat = [tmat; i+t_int];
    
end

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

%sgtitle('Multiparametric MPC')
sgtitle('Traditional MPC')