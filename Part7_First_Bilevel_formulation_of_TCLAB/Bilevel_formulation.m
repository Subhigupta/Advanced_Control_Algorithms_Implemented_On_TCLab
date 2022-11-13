clear all
clc

load problem.mat

% % Outer level formulation of objective function
% BPOP_problem.Q11=1;
% BPOP_problem.Q12=1;
% BPOP_problem.c11=zeros(1,1);
% BPOP_problem.c12=zeros(1,1);
% BPOP_problem.cc1=zeros(1,1);
% 
% % Outer level constraints
% BPOP_problem.A11=0;
% BPOP_problem.E11=0;

% Inner level formulation of objective function
BPOP_problem.Q22=problem.Q;
BPOP_problem.Q2b=problem.Ht;
BPOP_problem.c22=problem.c;

BPOP_problem.Q21=problem.Qt;
BPOP_problem.c21=problem.ct;
BPOP_problem.cc2=problem.cc;

% Inner level constraints
BPOP_problem.A22=problem.A;
BPOP_problem.b2=problem.b;
BPOP_problem.A21=-problem.F;
BPOP_problem.E21=zeros(26,1);
BPOP_problem.E22=zeros(26,1);

BPOP_problem.x1U=100;
BPOP_problem.x1L=0;
BPOP_Solution = BPOP(BPOP_problem);
