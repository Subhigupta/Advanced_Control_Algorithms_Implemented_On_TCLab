clear all
clc

problem.Q=[1 1 0;0 0.5 0;0 0 0.5];
problem.c=[1;1;0];
problem.Ht=[-3 0 0;0 1 0]';

problem.A=[2 1 -1;0 0 0;-1 -1 -1];
problem.b=[-2;0;0];
problem.F=[-1 2;1 1;0 0];

solution=mpQP(problem);

for i=1: Size(solution,2)
    H=0.5.*[1 0 -1;0 0 0;0 0 1];
    f=-4.*[0 1 0]-7.*

