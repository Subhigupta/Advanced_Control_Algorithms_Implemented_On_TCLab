# Model_Predictive_control_On_TCLab
The codes in this file demonstrate implementation of State-Space MPC, Explicit MPC and Bilevel MPC on an arduino device TCLab. 
The task performed by each of the code is described as follows:
1. Dynamic Modelling of TCLab
   - TCLab, a two-input two-output system, has been modeled as a dual heater system. The control inputs or manipulated variables are the heater outputs and the two controlled outputs are the temperature of the two heaters. The differential equations have been used to represent the system.
2. Linearising around a steady state point
   - The differential equations used to represent the system are a form of non-linear dynamic set of equations. The non-linear equations have been linearised around a steady-state point. The heaters were first kept on at 30% and were allowed to reach a steady state. A step change was given to capture the effect of input on output. The input-output data generated using step-response analysis was then fed to System Identification Toolbox, which as a result gave a linearised model, a state-space represention. This state-space representation is used as constraints in Linear Model Predictive Control formulation.
3. Case Study 1
     - In the first case study, the upper level objective aimed at minimizing a linear economic cost function while the lower level focused on tracking temperature of heater 1. 
4. Case Study 2
    - In the second case study, the upper level objective aimed at minimizing a linear economic cost function while the lower level focused on tracking temperature of heater 2. 
5. Case Study 3
    - In the third case study, the upper level objective focused tracking of temperature of heater 1 while the lower level focused on tracking temperature of heater 2.
In all the above case studies, offline solutions are generated for lower and upper level, offline solutions generation is a demonstration of Explicit MPC algorithm implementation.
The online procedure is used to calculate the parameter vector at each time instant and search of critical regions to get optimal control input.
