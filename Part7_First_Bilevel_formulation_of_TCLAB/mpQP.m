%==========================================================================
% File Name     : <mpQP.m>
% Usage         : [Solution,Time,Outer] = mpQP(problem,options)
% Description   : This function solves multi-parametric linear and
% quadratic programming problems based on the algorithm specified in
% OptionSet, unless specified otherwise in the optional input 'options'.
%--------------------------------------------------------------------------
% Author        : Richard Oberdieck, Nikolaos A. Diangelakis, 
%                 Efstratios N. Pistikopoulos
% Office        : Engineering Research Building, Texas A&M University, USA
% Mail          : paroc@tamu.edu
%--------------------------------------------------------------------------
% Last Revision | Author  | Description
%---------------+---------+------------------------------------------------
% 08-Mar-2016   | RO      | Initial Version
%---------------+---------+------------------------------------------------
% 26-Sep-2016   | RO      | Update to Version 2.0
%---------------+---------+------------------------------------------------
% 23-Feb-2017   | NAD     | Bug fixes
%---------------+---------+------------------------------------------------
% 18-Jun-2017   | NAD     | Bug fixes - ActiveSet index (reduced solution)
%---------------+---------+------------------------------------------------
% 12-Jul-2017   | NAD     | Bug fixes - ActiveSet index (full solution)
%==========================================================================

function [Solution,Time,Outer] = mpQP(problem,options)

% Define all the options
if nargin < 2
    options = OptionSet;
else
    Default = OptionSet;
    List_Fields = fieldnames(Default);
    for k = 1:length(List_Fields)
        field = List_Fields{k};
        if ~isfield(options,List_Fields{k}) || ...
                (isfield(options,field) && isempty(options.(field)))
            options.(field) = Default.(field);
        end
    end
end

% If mpQP is run, it is expected that they want to solve the mpQP part, so
% we remove the 'E'
if isfield(problem,'E')
    problem = rmfield(problem,'E');
    if isfield(problem,'Eeq')
        problem = rmfield(problem,'Eeq');
    end
end

[problem, ok] = ProblemSet(problem,options);

if ~ok
    Solution = [];
    Time = [];
    Outer = [];
    return
end

%% Go in the wrapper
switch lower(options.mpSolver)
    case 'mpt'
        [~, POP_Output, Tim] = POPviaMPT(problem);
        Solution = POP_Output.Solution;
        Outer = [];
        Time.Total = Tim;
        if ischar(Tim)
            Solution = 'TimedOut';
            Time.Total = Inf;
        end

        if strcmpi(options.SolutionStyle,'reduced')
            disp('The reduced solution is not available from MPT');
        end
        
    case 'geometrical'
        
        [Solution,Time,Outer] = Geometrical(problem,options);

    case 'combinatorial'
        
        [Solution,Time] = Combinatorial(problem,options);
        Outer = [];
        
    case 'graph'
        
        [Solution,Time] = ConnectedGraph(problem,options);
        Outer = [];
        
    case 'parallel_graph'
        % We switch off a specific warning
        warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
        [Solution,Time] = Parallel_ConnectedGraph(problem,options);
        Outer = [];
        
    case 'parallel_combinatorial'
        % We switch off a specific warning
        warning('off','MATLAB:mir_warning_maybe_uninitialized_temporary');
        [Solution,Time] = Parallel_Combinatorial(problem,options);
        Outer = [];
        
    otherwise
        
        disp('Unknown solver!');
        Solution = NaN;
        Time = NaN;
        Outer = NaN;
end
%% Here we change the indexing of the inequality constraints [for ActiveSet only]
if isfield(problem,'original') && ~isempty(Solution) && ~ischar(Solution)
    if isfield(Solution,'ActiveSet') % for full solution style
        for i = 1:size(Solution,2)
            Solution(i).ActiveSet = problem.original(1,Solution(i).ActiveSet(~isnan(Solution(i).ActiveSet)));
        end
    elseif ~strcmpi(options.mpSolver,'mpt') % for reduced solution style
        for i = 1:size(Solution,1)
            Solution(i,~isnan(Solution(i,:))) = ...
                problem.original(1,Solution(i,~isnan(Solution(i,:))));
        end
    end
end
end
