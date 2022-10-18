clear all
clc

load problem.mat

% Outer level formulation of objective function
BPOP_problem.Q11=1;
BPOP_problem.Q12=1;

% Outer level constraints
BPOP_problem.A11=0;
BPOP_problem.E11=0;
BPOP_problem.c21=0;
BPOP_problem.E21=0;

% Inner level formulation of objective function
BPOP_problem.Q22=problem.Q;
BPOP_problem.Q2b=problem.Ht;
BPOP_problem.Q21=problem.Qt;

% Inner level constraints
BPOP_problem.A22=problem.A;
BPOP_problem.b=problem.b;
BPOP_problem.A21=-problem.F;
BPOP_problem.E22=0;
BPOP_problem.cc2=0;
BPOP_problem.c22=0;
BPOP_problem.b2=0;

BPOP_Solution = BPOP(BPOP_problem);
