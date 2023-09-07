clear all
clc

%Inner level objective function
problem.Q22=[1 1 0;0 0.5 0;0 0 0.5];
problem.Q2b=[-3 0 0;0 1 0;0 0 0;0 0 0]';
problem.c22=[1;1;0];
problem.Q21=[0 0 0 0;0 0 0 0;0 0 0 0;0 0 0 0];
problem.c21=[0;0;0;0];
problem.cc2=0;

%Inner level constraints
problem.A22=[2 1 -1];
problem.b2=-2;
problem.A21=[1 -2 0 -1];
problem.E22=[];
problem.E21=[];

%Outer level objective function
problem.Q11=[0 0;0 0];
problem.Q1b=zeros(5,2)';
problem.c11=[-7;4];
problem.Q12=[1 0 -1 0 0;0 0 0 0 0;0 0 1 0 0;0 0 0 0 0;0 0 0 0 0];
problem.c12=[0;-4;0;0;0];
problem.cc1=0;

%Outer level constraints
problem.A11=[1 1];
problem.A12=[0 0 0 -1 0];
problem.E12=[];
problem.E11=[];
problem.b1=1;
problem.E11=[];

% Bounds
problem.x1U=[inf;inf];
problem.x1L=[0;0];
problem.x2U=[inf;inf;inf];
problem.x2L=[0;0;0];

Solution=BPOP(problem);