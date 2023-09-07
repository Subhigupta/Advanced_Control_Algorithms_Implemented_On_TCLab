clear all
clc

problem.Q=[1 0;0 1];
problem.Qt=[1 0;0 1];
problem.c=[40;40];
problem.ct=[-40;-40];
problem.cc=800;
problem.Ht=[-2 0;0 -2];

problem.A=[2 0;0 2;0 0;0 0;0 0;0 0;-1 0;0 -1;1 0;0 1];
problem.b=[-10;-10;0;0;50;50;10;10;20;20];
problem.F=[1 0;0 1;1 0;0 1;-1 0;0 -1;0 0;0 0;0 0;0 0];

Solution=mpQP(problem);
% Nine critical regions have been generated

for i=1:size(Solution,2)
    outer_obj=2*[1 0]+ 2*[0 1]-3.*Solution(i).Solution.X(1,1:end-1)-3.*Solution(i).Solution.X(2,1:end-1);
    outer_problem=outer_obj(:,1:end);
    A=Solution(i).CR.A;
    b=Solution(i).CR.b;
    [x,fval_outer_program]=cplexlp(outer_problem,A,b);
    final_sol(i).X=x;
    final_sol(i).Obj=fval_outer_program -3*Solution(i).Solution.X(1,end)-3*Solution(i).Solution.X(2,end)-60;
end