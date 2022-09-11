%    49.0000  314.7691  314.8556        % Centralized
%     49.0000  314.9783  315.0841       % Distributed


%%
clear 
clc
% close all

%%
%Load the matrix with the definition of the state space model

% load('ss1_Cen.mat');             %state space model step size 1
load('ss1.mat');  
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

%The MPC tuning Weights of MPC
%Just now just defined them as identity matrix
weighted_coefficient1_R1= 100; %0.001;
weighted_coefficient2_R1= 10; %0.01;
%mpc_crystal.R1 = 0.001*eye(size(model_crystal.B,2)); %Weight matrix for output moves (∆u)
mpc_crystal.R1 =blkdiag(weighted_coefficient1_R1,weighted_coefficient2_R1); %Weight matrix for output moves (∆u)
mpc_crystal.R=0.001*eye(size(model_crystal.B,2)); %Weight matrix for control inputs 
weighted_coefficient1_QR=100000;
weighted_coefficient2_QR=100000;
mpc_crystal.QR=blkdiag(weighted_coefficient1_QR,weighted_coefficient2_QR);
%mpc_crystal.QR = Weighted_Coefficients*eye(size(model_crystal.D,1));%Quadratic matrix for tracked output

%Input and Output Horizons
mpc_crystal.OH = 10; %10 %how many looks into the future for error minimzation
mpc_crystal.NC = 10;%10  % how many future control actions

%The bounds of the temperature states
mpc_crystal.Xmin = -1000000000*ones(size(model_crystal.A,1),1); 
mpc_crystal.Xmax = 1000000000*ones(size(model_crystal.A,1),1);

%The bounds of the inputs
mpc_crystal.Umin = [0 0]';   %lower bounds for u=Q1,Q2
mpc_crystal.Umax = [100 100]';    %upper bounds for u =Q1,Q2
mpc_crystal.Ymin=[0 0]';
mpc_crystal.Ymax=[100 100]'+273.15;

%Model-Plant mismatch (error term)
mpc_crystal.Ymismatch = [1;2];

problem = ss2qp_yalmip(mpc_crystal,model_crystal);