clear 
clc
%% Load the matrix with the definition of the state space model
%state space model step size 1

load('ss1.mat');  % used as a constraint in MPC formulation
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
model_crystal.D = ss1.C; % the matrix D in ss2qp is equivalent to matrix C in ss1

%% Matrix QR defined

weighted_coefficient1_QR = 5000;%10;%10;%1000;
weighted_coefficient2_QR = 5000;%10;%10;%1000
mpc_crystal.QR=blkdiag(weighted_coefficient1_QR,weighted_coefficient2_QR);
%mpc_crystal.QR = Weighted_Coefficients*eye(size(model_crystal.D,1));%Quadratic matrix for tracked output

%% Input and Output Horizons
% the Q matrix dimensions (NC*2) (equivalent to decision variables) depend on value of control horizon taken
% Since my manipulated variables are two then at each time step in the range of
% NC we will get two optimal manipulated variables

mpc_crystal.OH =2; %10 %how many looks into the future for error minimzation
mpc_crystal.NC =1;%10  % how many future control actions

%% Constraints

%The bounds of the temperature states
mpc_crystal.Xmin = -1000000*ones(size(model_crystal.A,1),1);%-10000000 
mpc_crystal.Xmax = 1000000*ones(size(model_crystal.A,1),1);%10000000

%The bounds of the inputs
mpc_crystal.Umin = [0 0]' - unominal';   %lower bounds for u=Q1,Q2
mpc_crystal.Umax = [100 100]'- unominal';    %upper bounds for u =Q1,Q2
mpc_crystal.Ymin=[0 0]'+273.15 - ynominal';
mpc_crystal.Ymax=[100 100]'+273.15- ynominal';

%The bounds on ∆u
% mpc_crystal.DUmax=[20 20]';
% mpc_crystal.DUmin=[-20 -20]';

%% Model-Plant mismatch (error term)
mpc_crystal.Ymismatch = [1;2];

%% Reformulating an optimal control into the quadratic programming problem using ss2qp
% ss2qp takes two arguments: state space model and definition of MPC control 
% (refer the word file on ss2qp or user manual on POP).

problem = ss2qp_yalmip(mpc_crystal,model_crystal);
% problem struct has information on variable and parameter vector.
problem.Q = problem.Q/2;

%% for mpMPC

%Type of solver for mp-LP/mp-QP. Default: ’Graph’.
%Uses the connected-graph algorithm.

%options.mpSolver='Graph';
%options.TimeMax =24*3600;
Solution=mpQP(problem);
%PlotSolution(Solution)
save Solution Solution
load('Solution.mat');