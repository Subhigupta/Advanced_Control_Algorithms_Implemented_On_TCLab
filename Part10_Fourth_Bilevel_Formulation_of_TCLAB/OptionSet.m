%==========================================================================
% File Name     : <OptionSet.m>
% Usage         : options = OptionSet
% Description   : This file was automatically created by POP.
% 
% This function contains the options set when the problem was exported. If
% the user would like to change these options, he or she is kindly referred
% to the user manual, where all the options are listed.
%--------------------------------------------------------------------------
% Author        : Richard Oberdieck, Nikolaos A. Diangelakis
%                 Efstratios N. Pistikopoulos
% Office        : Engineering Research Building, Texas A&M University, USA
% Mail          : paroc@tamu.edu
%--------------------------------------------------------------------------
% Generation date | Author  | Description
%-----------------+---------+----------------------------------------------
% 05-Jul-2020     | None    | Automatically generated version
%==========================================================================

function options = OptionSet
options.LPSolver = 'CPLEX';
options.QPSolver = 'CPLEX';
options.MISolver = 'CPLEX';
options.mpSolver = 'GRAPH';
options.SolutionStyle = 'full';
options.tolerance = sqrt(eps);
options.MinimalStep = 1000;
options.Step = 10;
options.TimeMax = 600;
options.BigM = 1e7;
options.Comparison = 'None';
options.DeltaSolver = 'Approximate';
options.rho_limit = 100;
options.Progress = 'text';
end