clear 
clc
%%
%Load the matrix with the definition of the state space model
%state space model step size 1

load('ss1.mat');  
load('unominal.mat');
load('ynominal.mat')
% %   Discrete-time identified state-space model:
% %     x(t+Ts) = A x(t) + B u(t) + K e(t)
% %        y(t) = C x(t) + D u(t) + e(t)

% % from ss2qp word document  
% % The usual linear discrete time state space model (LTI) equations:
% % 
% % x(t+1) 	= Ax(t) + Bu(t) + Cd(t)
% % y(t) 	= Dx(t) + Eu(t) + e

model_crystal.A = ss1.A;
model_crystal.B = ss1.B;
model_crystal.D = ss1.C;

%% 

%The MPC tuning Weights of MPC
weighted_coefficient1_R1 = 0.005;%0.1;%1; %100%0.001;
weighted_coefficient2_R1 = 0.005;%0.1;1;%0.1; %10 %0.01;
mpc_crystal.R1 =blkdiag(weighted_coefficient1_R1,weighted_coefficient2_R1); %Weight matrix for output moves (∆u)

%% 
weighted_coefficient_1 = 0.001/1000; %0.001
weighted_coefficient_2 = 0.001/1000; %0.001
mpc_crystal.R=blkdiag(weighted_coefficient_1,weighted_coefficient_2); %Weight matrix for control inputs 

%% 

weighted_coefficient1_QR = 10;%10;%1000;
weighted_coefficient2_QR = 10;%10;1000;
mpc_crystal.QR=blkdiag(weighted_coefficient1_QR,weighted_coefficient2_QR);

%% 

%Input and Output Horizons
mpc_crystal.OH = 2; %10 %how many looks into the future for error minimzation
mpc_crystal.NC = 2;%10  % how many future control actions
%% 

%The bounds of the temperature states
mpc_crystal.Xmin = -10000000*ones(size(model_crystal.A,1),1);%-10000000 
mpc_crystal.Xmax = 10000000*ones(size(model_crystal.A,1),1);%10000000

%The bounds of the inputs
mpc_crystal.Umin = [0 0]'-unominal';   %lower bounds for u=Q1,Q2
mpc_crystal.Umax = [100 100]' - unominal';    %upper bounds for u =Q1,Q2
mpc_crystal.Ymin=[0 0]'+273.15 - ynominal';
mpc_crystal.Ymax=[100 100]'+273.15 - ynominal';

%The bounds on ∆u
% mpc_crystal.DUmax=[20 20]';
% mpc_crystal.DUmin=[-20 -20]';

%% 

%Model-Plant mismatch (error term)
mpc_crystal.Ymismatch = [1;2];

problem = ss2qp_yalmip(mpc_crystal,model_crystal);
problem.Q = problem.Q/2;
%% 
%for mpMPC

%Type of solver for mp-LP/mp-QP. Default: ’Graph’.
%Uses the connected-graph algorithm.

% options.mpSolver='Graph';
% options.TimeMax =24*3600;
% Solution=mpQP(problem,options);
% % save Solution Solution
% load('Solution.mat');