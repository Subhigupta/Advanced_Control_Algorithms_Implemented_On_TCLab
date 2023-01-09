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

% % The below two equations define your model crystal formulation 
% % x(t+1) 	= Ax(t) + Bu(t) + Cd(t)
% % y(t) 	= Dx(t) + Eu(t) + e

model_crystal.A = ss1.A;
model_crystal.B = ss1.B(:,1); % in inner level u1 is kept as decision variable
model_crystal.C=ss1.B(:,2); %in inner level u2 is disturbance
model_crystal.D = ss1.C; % the matrix D in ss2qp is equivalent to matrix C in ss1

%% Matrix QR defined : Weights for temperature set point tracking

weighted_coefficient1_QR = 1500000;%10;%10;%1000;
weighted_coefficient2_QR = 0;%10;%10;%1000
mpc_crystal.QR=blkdiag(weighted_coefficient1_QR,weighted_coefficient2_QR);
%Quadratic matrix for tracked output
mpc_crystal.R=50;
mpc_crystal.R1=50000;
%% Input and Output Horizons
% the Q matrix dimensions (NC*2) (equivalent to decision variables) depend on value of control horizon taken
% Since my manipulated variables are two then at each time step in the range of
% NC we will get two optimal manipulated variables

mpc_crystal.OH = 2; %10 %how many looks into the future for error minimzation
mpc_crystal.NC = 1;%10  % how many future control actions

%% Constraints

%The bounds of the temperature states
mpc_crystal.Xmin = -1000000*ones(size(model_crystal.A,1),1);%-10000000 
mpc_crystal.Xmax = 1000000*ones(size(model_crystal.A,1),1);%10000000

%The bounds of the inputs
%mpc_crystal.Umin = 0 - unominal';
mpc_crystal.Umin = 0 - unominal';%lower bounds for u=Q1,Q2
mpc_crystal.Umax = 100- unominal';    %upper bounds for u =Q1,Q2

% mpc_crystal.Dmin = 0- unominal';   %lower bounds for disturbance
mpc_crystal.Dmin = 0- unominal';
mpc_crystal.Dmax = 100- unominal';    %upper bounds for disturbance

%So range of manipulated variable is -30 < u < 70

mpc_crystal.Ymin=[0 0]'+273.15 - ynominal';
mpc_crystal.Ymax=[100 100]'+273.15- ynominal';
% So range of controlled variable is -41.29<y1<58.70
% So range of controlled variable is -37.99<y2<62.00376

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

save Solution Solution
load('Solution.mat');

%% Solving outer level Second iteration
Solution_Inner_Level=[];
Solution_Outer_level=[];

for i=1:size(Solution,2)
    %Find dependent solution of inner problem in term of independent parameters
    Solution_Inner_Level=[Solution_Inner_Level;Solution(i).Solution.X];
    n=size(Solution(i).Solution.X,2);
    m=size(Solution(i).Solution.X,1);
    dependent_variable=(Solution_Inner_Level(i,1:n-1));
    
    %Find new objective function for outer problem
    u1=dependent_variable(:,5);
    u2=1;
    obj_new_outer=u1+u2;
    
    %Formulation of objective function of outer problem
    outer_problem.c=obj_new_outer;
    outer_problem.ct=[dependent_variable(1:4), dependent_variable(6:end)]';

    %Constraints of outer problem
    outer_problem.A=Solution(i).CR.A(:,5);
    outer_problem.b=Solution(i).CR.b;
    outer_problem.F=-Solution(i).CR.A(:,[1:4,6:end]);
 
    %Solving outer problem multiparametrically
    sol=mpQP(outer_problem);
    FinalSolution(i).ans=sol;
    for j=1:size(FinalSolution(i).ans,2)
        Solution_Outer_level=[Solution_Outer_level;FinalSolution(i).ans(j).Solution.X];  
    end
end
%% FinalSolution separated into critical regions and solution vector 

t=1;
for i=1:size(FinalSolution,2)
    for j=1:size(FinalSolution(i).ans,2)
        Final(t).CR=FinalSolution(i).ans(j).CR;
        Final(t).Solution=FinalSolution(i).ans(j).Solution;
        t=t+1;
    end
end