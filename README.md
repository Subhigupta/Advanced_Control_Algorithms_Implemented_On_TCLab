# Model_Predictive_control_On_TCLab
The codes in this file demonstrate implementation of State-Space MPC, Explicit MPC and Bilevel MPC on an arduino device TCLab. 
The task performed by each of the code is described as follows:
1. Dynamic Modelling of TCLab
   - TCLab is a two-input two-output system. It has been modeled as a dual heater system. The control inputs or manipulated variables are the heater outputs and the temperature of the two heaters are the two controlled outputs. The differential equations have been used to represent the system.
3. Linearising around a steady state point
   - The differential equations used to represent the system are a form of non-linear dynamic set of equations. The non-linear equations have been linearised around a steady-state point. The heaters were kept on at 30% and were allowed to reach a steady state. A step change was given to capture the effect of input on output. The input-output data generated used then fed to system identification toolbox, which as a result gave a linearised model, a state-space represention.
5. Case Study 1
     - In the first case study, the upper level objective is minimization of a linear economic cost function while the lower level focused on tracking temperature of heater 1. 
6. Case Study 2
    - In the second case study, the upper level objective is minimization of a linear economic cost function while the lower level focused on tracking temperature of heater 2. 
7. Case Study 3
    - In the third case study, the upper level objective focuses on tracking of temperature of heater 1 while the lower level focused on tracking temperature of heater 2. 
