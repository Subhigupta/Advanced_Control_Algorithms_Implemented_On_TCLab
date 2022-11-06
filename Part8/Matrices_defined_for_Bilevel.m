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
model_crystal.B = ss1.B(:,2);
model_crystal.C = ss1.B(:,1); % keeping u1 as measured disturbance for inner level decision
% variable being u2
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
mpc_crystal.Umin = 0;  %lower bounds for u=Q1,Q2
mpc_crystal.Umax = 100;    %upper bounds for u =Q1,Q2

% mpc_crystal.Dmin = 0;   %lower bounds for u=Q1,Q2
% mpc_crystal.Dmax = 100;    %upper bounds for u =Q1,Q2

%So range of manipulated variable is -30 < u < 70
mpc_crystal.Ymin=[0 0]'+273.15;
mpc_crystal.Ymax=[100 100]'+273.15;
% So range of controlled variable is 41.29<y1<58.70
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
% tfixed = NaN*ones(size(Solution(1).CR.A,2),1);
% tfixed(1)=10000;
% tfixed(2)=10000;
% tfixed(3)=10000;
% tfixed(4)=10000;
% 
% tfixed(6)=40;
% 
% tfixed(8)=50;
% tfixed(9)=50;

%option = {'CR','OBJ','all'};
%PlotSolution(Solution,tfixed)%,tfixed)
save Solution Solution
load('Solution.mat');
%% Solving outer level first iteration

% for i=1:11
%     %Find dependent solution of inner problem in term of independent parameters
%     matrix_for_dependent_variable=Solution(i).Solution.X;
%     n=size(Solution(i).Solution.X,2);
%     m=size(Solution(i).Solution.X,1);
%     dependent_variable=(matrix_for_dependent_variable(1:m,1:n));
%     %dependent_variable=((matrix_for_dependent_variable(1:3,1:2))*[x1;x2])+ matrix_for_dependent_variable(1:3,3);
%     %Find new objective function for outer problem
%     %u2=dependent_variable(1,1:n);
%     weights=[10e-05 10e-05 10e-05 10e-05 1 10e-05 10e-05 10e-05 10e-05 10e-05];
%     u2=weights.*dependent_variable(1,1:n);
%     u1=zeros(1,size(u2,2));
%     u1(5)=1;
%     Aineq=Solution(i).CR.A;
%     bineq=Solution(i).CR.b;
%     obj_new_outer=u1+u2;
%    [x,fval_outer_program]=cplexlp(obj_new_outer(1,1:n-1),Aineq,bineq);
%    final_sol(i).X=x;
%    final_sol(i).Obj=fval_outer_program + obj_new_outer(1,10);
% end
%% Solving outer level Second iteration

for i=1:11
    %Find dependent solution of inner problem in term of independent parameters
    matrix_for_dependent_variable=Solution(i).Solution.X;
    n=size(Solution(i).Solution.X,2);
    m=size(Solution(i).Solution.X,1);
    dependent_variable=(matrix_for_dependent_variable(1:m,1:n-1));
    
    %Find new objective function for outer problem
    u2=zeros(m,n-1);
    u2(:,5)=dependent_variable(:,5);
    u1=zeros(1,size(u2,2));
    u1(:,5)=1;
    obj_new_outer=u1+u2;
    
    %Formulation of objective function of outer problem
    outer_problem.c=obj_new_outer';
    outer_problem.ct=[dependent_variable(1:4),0, dependent_variable(6:n-1)]';

    %Constraints of outer problem
    A_matrix=zeros(size(Solution(i).CR.A,1),size(Solution(i).CR.A,2));
    A_matrix(:,5)=Solution(i).CR.A(:,5);
    F_matrix=zeros(size(Solution(i).CR.A,1),size(Solution(i).CR.A,2));
    F_matrix(:,1:4)=Solution(i).CR.A(:,1:4);
    F_matrix(:,6:end)=Solution(i).CR.A(:,6:end);
    outer_problem.A=A_matrix;
    outer_problem.b=Solution(i).CR.b;
    outer_problem.F=-F_matrix;
    
    %Solving outer problem multiparametrically
    sol=mpQP(outer_problem);
    FinalSolution(i).ans=sol;
end