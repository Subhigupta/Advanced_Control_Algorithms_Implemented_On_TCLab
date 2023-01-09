clear 
clc

%% Load the matrix with the definition of the state space model
%state space model step size 1

load('ss1.mat'); 

model_crystal_inner.A = ss1.A;
model_crystal_inner.B = ss1.B(:,2); % in inner level u2 is kept as decision variable
model_crystal_inner.C=ss1.B(:,1); %in inner level u1 is disturbance
model_crystal_inner.D = ss1.C; % the matrix D in ss2qp is equivalent to matrix C in ss1

%% Inner Level Formulation (Tracking T2)

weighted_coefficient1_QR = 0;
weighted_coefficient2_QR = 20000000; 
mpc_crystal_inner.QR=blkdiag(weighted_coefficient1_QR,weighted_coefficient2_QR);
mpc_crystal_inner.R = 200;

mpc_crystal_inner.OH = 2;
mpc_crystal_inner.NC = 1;

%The bounds of the temperature states
mpc_crystal_inner.Xmin = -1000*ones(size(model_crystal_inner.A,1),1);%-10000000 
mpc_crystal_inner.Xmax = 1000*ones(size(model_crystal_inner.A,1),1);%10000000

%The bounds of the inputs
mpc_crystal_inner.Umin = 0 - unominal';
%mpc_crystal.Umin = 28 - unominal';%lower bounds for u=Q1,Q2
mpc_crystal_inner.Umax = 100- unominal';    %upper bounds for u =Q1,Q2

mpc_crystal_inner.Dmin = 0- unominal';   %lower bounds for disturbance
mpc_crystal_inner.Dmax = 100- unominal';    %upper bounds for disturbance
%So range of manipulated variable is -30 < u < 70

mpc_crystal_inner.Ymin=[0 0]'+273.15 - ynominal';
mpc_crystal_inner.Ymax=[100 100]'+273.15- ynominal';
% So range of controlled variable is -41.29<y1<58.70
% So range of controlled variable is -37.99<y2<62.00376

mpc_crystal_inner.Ymismatch = [1;2];

problem_inner = ss2qp_yalmip(mpc_crystal_inner,model_crystal_inner);
problem_inner.Q = problem_inner.Q/2;

% Solving inner level multiparamet
Solution_inner=mpQP(problem_inner); 

save Solution_inner Solution_inner
load('Solution_inner.mat');

%% Outer Level Formulation (Tracking T1)

model_crystal_outer.A = ss1.A;
model_crystal_outer.B = ss1.B(:,1); % in outer level u1 is kept as decision variable
model_crystal_outer.C=ss1.B(:,2); %in outer level u2 is disturbance
model_crystal_outer.D = ss1.C; % the matrix D in ss2qp is equivalent to matrix C in ss1

weighted_coefficient1_QR = 20000000;%1000;
weighted_coefficient2_QR = 0;%10;%10;%1000
mpc_crystal_outer.QR=blkdiag(weighted_coefficient1_QR,weighted_coefficient2_QR);
mpc_crystal_outer.R = 900;%10;

mpc_crystal_outer.OH = 2;
mpc_crystal_outer.NC = 1;

%The bounds of the temperature states
mpc_crystal_outer.Xmin = -1000*ones(size(model_crystal_outer.A,1),1);%-10000000 
mpc_crystal_outer.Xmax = 1000*ones(size(model_crystal_outer.A,1),1);%10000000

%The bounds of the inputs
mpc_crystal_outer.Umin = 0 - unominal';
%mpc_crystal.Umin = 28 - unominal';%lower bounds for u=Q1,Q2
mpc_crystal_outer.Umax = 100- unominal';    %upper bounds for u =Q1,Q2

mpc_crystal_outer.Dmin = 0- unominal';   %lower bounds for disturbance
mpc_crystal_outer.Dmax = 100- unominal';    %upper bounds for disturbance
%So range of manipulated variable is -30 < u < 70

mpc_crystal_outer.Ymin=[0 0]'+273.15 - ynominal';
mpc_crystal_outer.Ymax=[100 100]'+273.15- ynominal';

mpc_crystal_outer.Ymismatch = [1;2];

problem_outer = ss2qp_yalmip(mpc_crystal_outer,model_crystal_outer);
problem_outer.Q = problem_outer.Q/2;

%% Reformulating outer level quadratic problem after substituting solution and critical region of inner level

for i=1:size(Solution_inner,2)
    problem_outer_BLPP.Ht=problem_outer.Ht(:,[1:4,6:end])+ (2*problem_outer.Qt(5,[1:4,6:end]))*(Solution_inner(i).Solution.X(5));
    duplicate_Qt=problem_outer.Qt;
    duplicate_Qt(:,5)=[];
    duplicate_Qt(5,:)=[];
    problem_outer_BLPP.Qt=duplicate_Qt;
    problem_outer_BLPP.Q=problem_outer.Q + (((Solution_inner(i).Solution.X(5))^2)*problem_outer.Qt(5,5)) + (problem_outer.Ht(:,5)*Solution_inner(1).Solution.X(5));
    problem_outer_BLPP.A=Solution_inner(i).CR.A(:,5);
    problem_outer_BLPP.b=Solution_inner(i).CR.b;
    problem_outer_BLPP.F=-Solution_inner(i).CR.A(:,[1:4,6:end]);
   
    sol=mpQP(problem_outer_BLPP);
    FinalSolution(i).ans=sol;
end
%% FinalSolution separated into critical regions and solution vector 

t=1;
for i=1:size(FinalSolution,2)
    for j=1:size(FinalSolution(i).ans,2)
        Solution_outer(t).CR=FinalSolution(i).ans(j).CR;
        Solution_outer(t).Solution=FinalSolution(i).ans(j).Solution;
        t=t+1;
    end
end

save Solution_outer Solution_outer
load('Solution_outer.mat');