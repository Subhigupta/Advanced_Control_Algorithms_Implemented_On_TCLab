clear
clc

%% Problem Formulation

problem.Q11 = blkdiag(4, 1, 0);
problem.Q1b = zeros(3);
problem.c11 = [0; 2; 6];
problem.Q12 = zeros(3);
problem.c12 = [1; 5; 0];
problem.cc1 = [0];

problem.Q22 = blkdiag(4, 1, 5);
problem.Q2b = [0, 0, 0; 0, 1, 0; 0, -1, 0];
problem.c22 = [-5; -15; -16];
problem.Q21 = zeros(3);
problem.c21 = zeros(3, 1);
problem.cc2 = [0];

problem.A11 = [0, 0];
problem.E11 = [-1];
problem.A12 = [0];
problem.E12 = [-1, 1];
problem.b1 = [-1];
problem.A21 = [[6.4, 7.2; -8, -4.9; 3.3, 4.1]; zeros(6, 2)];
problem.E21 = [0; 0; 0.5; -1; 1; zeros(4, 1)];
problem.A22 = [2.5; -3.2; 0.02; zeros(6, 1)];
problem.E22 = [zeros(2); [4, 4.5]; zeros(2); -eye(2); eye(2)];
problem.b2 = [11.5; 5; 1; 0; 1; 0; 0; 1; 1];
problem.x1L = [-10; -10];
problem.x1U = [10; 10];


%% User settings

options.LPSolver = 'CPLEX';
options.QPSolver = 'CPLEX';
options.MISolver = 'CPLEX';
options.mpSolver = 'Combinatorial';
options.SolutionStyle = 'full';
options.Comparison = 'Exact';
ModifyOptionSet(options);


%% Solve the problem

[out_BPOP, out_POP] = BPOP(problem)