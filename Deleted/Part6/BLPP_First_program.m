clear all
clc

%syms x1 x2 y1 y2 y3
% obj_outer=-8*x1-4*x2+4*y1-40*y2+4*y3;
% obj_inner=x1+2*x2+y1+y2+2*y3;

%Make a struct of your inner problem
BLPP.c=[1;1;2];
BLPP.ct=[1;2];
BLPP.A=[-1 1 1;-1 2 -0.5;2 -1 -0.5;-1 0 0;0 -1 0;0 0 -1;0 0 0;0 0 0];
BLPP.b=[1;1;1;0;0;0;0;0];
BLPP.F=[0 0;-2 0;0 -2;0 0;0 0;0 0;1 0;0 1];

%Solve the inner problem using multi-parametric technique.
sol=mpQP(BLPP);

%Plot the parametric solution 
PlotSolution(sol)

for i=1:5
    %Find dependent solution of inner problem in term of independent parameters
    matrix_for_dependent_variable=sol(i).Solution.X;
    dependent_variable=(matrix_for_dependent_variable(1:3,1:3));
    %dependent_variable=((matrix_for_dependent_variable(1:3,1:2))*[x1;x2])+ matrix_for_dependent_variable(1:3,3);
    %Find new objective function for outer problem
    y1=dependent_variable(1,1:3);
    y2=dependent_variable(2,1:3);
    y3=dependent_variable(3,1:3);
    x1=[1 0 0];
    x2=[0 1 0];
    obj_new_outer=-8*x1-4*x2+4*y1-40*y2+4*y3;
    A=sol(i).CR.A;
    b=sol(i).CR.b;
   [x,fval_outer_program]=cplexlp(obj_new_outer(1,1:2),A,b);
   final_sol(i).X=x;
   final_sol(i).Obj=fval_outer_program + obj_new_outer(1,3);
end
