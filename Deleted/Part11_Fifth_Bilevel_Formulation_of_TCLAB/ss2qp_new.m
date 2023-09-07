%==========================================================================
% File Name     : <ss2qp.m>
% Usage         : [mpv, online] = ss2qp(mpc, ss)
% Description   : This function converts a state-space MPC controller into
% a mp-QP. It returns the mp-QP equivalent and the online version of the
% MPC controller. The online version can be used as a quick "parametric"
% MPC and also has information on variable and parameter vector.
%
%
% order of variable vector (dependent stuff on the left)
% each subvector ordered according to state-space matrices;
%                                         earlier time ones start from left
%
% The error terms are (Y - Yref) as they appear in the objective, they are
%         NOT aligned with Yi but however they appear in objective function
%
% NEW! "error" terms for U, offering reference points for controls
%                                                               Eu = U-Uref
% the last optional variable (scalar) is introduced when we have constraint
%                                             slackening (p.99 Maciejowski)
%
% variables: [Du0...DuN-1(deltas) y1...yN-1(outputs) x1...xN(states)
%    e1...eN-1(errors) Eu0...EuN-1(control error) u0...uN-1(controls) soft]
%
% Du, Y and Error variables are optional depending on the configuration
%
% naming convention: X=state, U=control, D=disturbance, Y=output
%--------------------------------------------------------------------------
% Authors       : Richard Oberdieck, Nikolaos Diangelakis, Christos Panos,
%                 Nikolaos Bozinis, Efstratios N. Pistikopoulos
% Office        : Engineering Research Building, Texas A&M University, USA
% Mail          : paroc@tamu.edu
%--------------------------------------------------------------------------
% Last Revision | Author  | Description
%---------------+---------+------------------------------------------------
% 31-Oct-2015   | RO      | Initial Version
%==========================================================================

function [mpv, online] = ss2qp(mpc, ss)

%% Initialization
% We define this like this for now. Maybe later we integrate it with the
% POP options completely
options = OptionSet;
ZERO = options.tolerance;

nuU = size(ss.B,2);
nuX = size(ss.A,2);
if ~isfield(ss, 'C')
    ss.C = zeros(nuX,0);
end
nuD = size(ss.C,2);
if ~isfield(ss, 'D')
    nuY = 0;
else
    nuY = size(ss.D,1);
end
nuE = 0; % "error" variables (Y-Yref)

if mpc.OH < mpc.NC || mpc.OH < 1 || mpc.NC <= 0
    error('Wrong horizon definition, please modify');
end

if nuU == 0 || nuX == 0
    error('No control and/or state variables detected');
end

if (isfield(mpc, 'R1') && isempty(mpc.R1)==0 && norm(mpc.R1, inf) > eps) || ...
        isfield(mpc, 'Dumax') || isfield(mpc, 'Dumin')
    nudU = nuU;
    % @@@ what if only a subset of controls are step-blocked?
else
    nudU = 0;
end

if nuY
    % if we have no E matrix, make sure it is zeroed out
    if ~isfield(ss, 'E') || isempty(ss.E)
        ss.E = zeros(nuY, nuU);
    end
end

% @@@ treat all mpc vectors to homogenous form (over horizon)


%% ----------- EQUATIONS: MPC LOOKAHEAD DEFINITIONS -----------
% all equations are implicit LHS = 0 (no explicit constant term)
% the column order is first variables, then group#1 of thetas
if isfield(mpc, 'bFixedDisturbance') && mpc.bFixedDisturbance
    nuT1Dis = nuD; % disturbance assumed fixed at D0 levels
else
    nuT1Dis = nuD*mpc.OH; % separate parameters for each prediction interval
end
if isfield(mpc, 'Ymismatch')
    if any(mpc.Ymismatch < 1) || any(mpc.Ymismatch > nuY) || ...
            length(mpc.Ymismatch) > length(unique(mpc.Ymismatch))
        error('output mismatch index Ymismatch error');
    end
    nT1Ymis = length(mpc.Ymismatch);
    % calculating the starting mismatch may require DU formulation (NO yer plonker! we already have u0)
    % if nT1Ymis & 0==nudU & norm(ss.E, 'inf') > eps, nudU = nuU; end
else
    nT1Ymis = 0;
end

nuT1 = nuX + nudU + nuT1Dis + nT1Ymis;
eqs = [];
namesVars = []; % all inclusive variable vector (dependents + independents)

% start with DU equations, if any [Du = u(t) - u(t-1)]; remember that time goes up left to right!
% yalmip me @rse <g>
if nudU
    lepad = zeros(nudU, nudU);
    stub = [eye(nudU) zeros(nudU, nudU*(mpc.OH-2) + nuY*(mpc.OH-1) + nuX*mpc.OH) eye(nudU) -eye(nudU)];
    ripad = zeros(nudU, nuU*(mpc.OH-2) + nuT1);
    % the first one incorporates U-1 parameter (past control action)
    eqs = [eqs; eye(nudU) zeros(nudU, nudU*(mpc.OH-1) + nuY*(mpc.OH-1) + nuX*mpc.OH) -eye(nuU) zeros(nudU, nuU*(mpc.OH-1) + nuX) eye(nuU) zeros(nudU, nuT1Dis+nT1Ymis)];
    for i=1:mpc.OH-1
        eqs = [eqs; lepad stub ripad];
        lepad = [lepad zeros(nudU, nudU)];
        ripad(:, 1:nudU) = [];
    end
    
    aux = 1;
    for i=1:mpc.OH
        for j=1:nuU
            namesVars{aux} = sprintf('Du%d(t+%d)', j, i-1);
            aux = aux+1;
        end
    end
end

% output definitions, from N=1 (ie not for the first interval)
% y = Dx + Eu + error
lepad = zeros(nuY, nudU*mpc.OH);
ripad = zeros(nuY, nuT1);
for i=1:mpc.OH-1
    if ~nuY
        break;
    end
    % note below how we jump to u1 (skip u0)
    stub = [eye(nuY) zeros(nuY, nuY*(mpc.OH-i-1)) zeros(nuY, nuX*(i-1)) -ss.D zeros(nuY, nuX*(mpc.OH-i)+nuU) zeros(nuY, nuU*(i-1)) -ss.E zeros(nuY, nuU*(mpc.OH-1-i))];
    eqs = [eqs; lepad stub ripad];
    lepad = [lepad zeros(nuY, nuY)];
    
    for j=1:nuY
        namesVars{length(namesVars)+1} = sprintf('y%d(t+%d)', j, i);
    end
    
    if nT1Ymis
        % error term is Yreal0 - Yo = Yr0 - (Dx0+Eu0) and that has to be added to all future Ys
        % (it is assumed constant for all future predictions)
        aux = (nudU+nuY+nuX)*mpc.OH-nuY + 1; % u0 column
        tmp = aux + mpc.OH*nuU; % first parameter
        idx = tmp + nuT1 - nT1Ymis;
        for j=mpc.Ymismatch
            row = size(eqs, 1) - nuY + j;
            
            eqs(row, aux:aux+nuU-1) = ss.E(j,:);
            eqs(row, tmp:tmp+nuX-1) = ss.D(j,:);
            eqs(row, idx) = -1;
            idx = idx + 1;
        end
    end
end

% state space dynamic equations: x+1 = Ax + Bu + Cd
% first one differs a lot
lepad = zeros(nuX, nudU*mpc.OH+nuY*(mpc.OH-1));
ripad = zeros(nuX, nuD*(mpc.OH-1) + nT1Ymis);
% local equation array assuming fully parametrized disturbance
eqtmp = [lepad eye(nuX) zeros(nuX, nuX*(mpc.OH-1)) -ss.B zeros(nuX, nuU*(mpc.OH-1)) -ss.A zeros(nuX,nudU) -ss.C ripad];
for i=1:mpc.OH-1
    stub = [-ss.A eye(nuX) zeros(nuX, nuX*(mpc.OH-i-1) + nuU*i) -ss.B zeros(nuX, nuU*(mpc.OH-i-1)+nuX+nudU+nuD) zeros(nuX, nuD*(i-1)) -ss.C];
    eqtmp = [eqtmp; lepad stub zeros(nuX, nuD*(mpc.OH-i-1) + nT1Ymis)];
    lepad = [lepad zeros(nuX, nuX)];
    for j=1:nuX
        namesVars{length(namesVars)+1} = sprintf('x%d(t+%d)', j, i);
    end
end
i=mpc.OH;
for j=1:nuX
    namesVars{length(namesVars)+1} = sprintf('x%d(t+%d)', j, i);
end

% if we have fixed disturbance policy, lose the extra columns
if nuT1Dis == nuD & nuD
    aux = (nudU+nuY+nuX+nuU)*mpc.OH -nuY +nuX + nudU + nuD + 1;
    eqtmp(:, aux:aux+(mpc.OH-1)*nuD-1) = [];
    aux = aux - nuD;
    
    % all future Xs depend on the starting disturbance D0
    eqtmp(:, aux:aux+nuD-1) = repmat(-ss.C, mpc.OH, 1);
end
eqs = [eqs; eqtmp];

for i=1:mpc.OH
    for j=1:nuU
        namesVars{length(namesVars)+1} = sprintf('u%d(t+%d)', j, i-1);
    end
end

% if control horizon is smaller, fix last few controls
aux = mpc.OH - mpc.NC;
if aux > 0
    % @@@ add u=Kx equation if we have stabilizing feedback
    lepad = zeros(nuU, (nudU+nuY+nuX+nuU)*mpc.OH -nuY -(aux+1)*nuU);
    stub = [eye(nuU) -eye(nuU)];
    ripad = zeros(nuU, nuU*(aux-1) + nuT1);
    for i=1:aux
        eqs = [eqs; lepad stub ripad];
        if i < aux
            lepad = [lepad zeros(nuU, nuU)];
            ripad(:, 1:nuU) = [];
        end
    end
end

% to minimize on control actions, we could condensate successive intervals
if isfield(mpc, 'uzip')
    if isfield(mpc, 'u_bundle')
        error('It does not make sense to set ''uzip'' and ''u_bundle''. Please change.');
    elseif mpc.uzip <= 1 || mpc.uzip > mpc.NC
        error('The definition of uzip is incorrect.');
    end
    
    % remember that index 1 corresponds to U0
    % this is short-cut definition
    mpc.u_bundle = [];
    for i=1:floor(mpc.NC/mpc.uzip)
        aux = (i-1)*mpc.uzip + 1;
        mpc.u_bundle{i} = aux:aux+mpc.uzip-1;
    end
    aux = aux + mpc.uzip;
    if aux <= mpc.NC
        mpc.u_bundle{i+1} = aux:mpc.NC; % odds and sods
    end
end
bUBundle = isfield(mpc, 'u_bundle');
if bUBundle
    
    aux = [];
    for j=1:length(mpc.u_bundle)
        aux = [aux mpc.u_bundle{j}];
    end
    % unsorted ranges are possible, albeit not too reasonable from a physical point of view!
    if length(aux) > length(unique(aux)) || sum(aux) ~= 0.5*mpc.NC*(mpc.NC+1)
        keyboard
        error('u bundling must cover all control horizon range in order');
    end
    % at the limit length(u_bundle)==NC where we don't really group anything
    
    % generate equations for each group
    % we start counting from the closest U in time towards HC
    eqsp = [];
    base = eye(nuU);
    for i=1:length(mpc.u_bundle)
        idx = sort(mpc.u_bundle{i});
        if length(idx) > 1
            ind = idx(1);
            lepad = zeros(nuU,(ind-1)*nuU);
            for k=2:length(idx)
                % add equation block U1 == U2
                i2 = idx(k);
                tmp = i2-ind-1; % how far apart these are, >=0! (since idx(1) is smaller)
                tmp = [base zeros(nuU,tmp*nuU) -base];
                eqsp = [eqsp; lepad tmp zeros(nuU, nuU*mpc.NC-size(lepad,2)-size(tmp,2))];
            end
        end
    end
    
    bUBundle = size(eqsp,1); % couples as extra dependent vars
    if length(eqsp)
        aux = (nudU+nuY+nuX)*mpc.OH - nuY;
        eqs = [eqs; zeros(bUBundle, aux) eqsp zeros(bUBundle, size(eqs,2)-aux-size(eqsp,2))];
    end
end

namesT1 = [];
aux = 1;
for j=1:nuX
    namesT1{aux} = sprintf('x%d(t)', j);
    aux = aux+1;
end
if nudU
    for j=1:nuU
        namesT1{aux} = sprintf('u%d(t-1)', j);
        aux = aux+1;
    end
end
if nuT1Dis == nuD
    idx = 1;
else
    idx = 1:mpc.OH;
end
for i=idx
    for j=1:nuD
        namesT1{aux} = sprintf('d%d(t+%d)', j, i-1);
        aux = aux+1;
    end
end
if nT1Ymis
    for i = mpc.Ymismatch
        namesT1{aux} = sprintf('y%d(t)_real', i);
        aux = aux+1;
    end
end

% how many dependent variables we will be eliminating
% a single extra variable will have to be added for constraint slackening (totally free, only appears in constraints)
ndep = (nudU+nuY+nuX)*mpc.OH -nuY + bUBundle + (mpc.OH - mpc.NC)*nuU;

%% ----------- CONSTRAINTS: BOUNDS ON VARIABLES -----------
% work in subsystems of T3: states, outputs etc
% all constraints in the form: X - T3 < rhs
% Additions by Richard Oberdieck for arbitrary constraint sets
conXM = [];
namesT3 = [];
nT3Xmax = 0;
if isfield(mpc, 'Xmax') && isfield(mpc, 'Xmin')
    aux = find(mpc.Xmax < mpc.Xmin);
    if length(aux) && isfinite(mpc.Xmax(aux)) && isfinite(mpc.Xmin(aux))
        error('The bounds on X are infeasible')
    end
end
if isfield(mpc, 'Xnocon'), nocon = mpc.Xnocon;
else nocon = [];
end

if isfield(mpc,'XConstA')
    if size(mpc.XConstA,2) == nuX
        XConstA = mpc.XConstA;
        for j = 2:mpc.OH
            mpc.XConstA = blkdiag(mpc.XConstA,XConstA);
        end
        mpc.XConstb = repmat(mpc.XConstb,mpc.OH,1);
    end    
    
    len = size(mpc.XConstA,1);
    sp1 = lower(len/2);
    conXM = [mpc.XConstA(1:sp1,:), mpc.XConstb(1:sp1)];
    conXL = [mpc.XConstA(sp1+1:end,:), mpc.XConstb(sp1+1:end)];
    nT3Xmax = 0;
    nT3Xmin = 0;
else
    
    if isfield(mpc, 'Xmax')
        if length(mpc.Xmax) ~= nuX, error('Xmax dimension foul'); end
        
        [conXM remain] = addbounds(mpc.Xmax, 1, mpc.OH, nocon);
        
        for j = 1:nuX
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('xMax%d', j);
            end
        end
        nT3Xmax = length(remain);
    end

    conXL = [];
    nT3Xmin = 0;
    if isfield(mpc, 'Xmin')
        if length(mpc.Xmin) ~= nuX, error('Xmin dimension foul'); end
        
        [conXL remain] = addbounds(mpc.Xmin, 0, mpc.OH, nocon);
        
        for j = 1:nuX
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('xMin%d', j);
            end
        end
        nT3Xmin = length(remain);
    end
end
conUM = [];
nT3Umax = 0;
if isfield(mpc, 'Unocon'), nocon = mpc.Unocon;
else nocon = [];
end
if isfield(mpc, 'Umax') & isfield(mpc, 'Umin')
    aux = find(mpc.Umax < mpc.Umin);
    if length(aux) & isfinite(mpc.Umax(aux)) & isfinite(mpc.Umin(aux))
        error('bounds on U are infeasible')
    end
end

if isfield(mpc,'UConstA')
    if size(mpc.UConstA,2) == nuU
        UConstA = mpc.UConstA;
        for j = 2:mpc.OH
            mpc.UConstA = blkdiag(mpc.UConstA,UConstA);
        end
        mpc.UConstb = repmat(mpc.UConstb,mpc.OH,1);
    end
    
    len = size(mpc.UConstA,1);
    sp1 = lower(len/2);
    conUM = [mpc.UConstA(1:sp1,:), mpc.UConstb(1:sp1)];
    conUL = [mpc.UConstA(sp1+1:end,:), mpc.UConstb(sp1+1:end)];
    nT3Umax = 0;
    nT3Umin = 0;
else
    if isfield(mpc, 'Umax')
        if length(mpc.Umax) ~= nuU, error('Umax dimension foul'); end
        
        % note: constraints on U go only up to the control horizon (no OH)
        % @@@ what about eliminating those for block-policy controls?
        [conUM remain] = addbounds(mpc.Umax, 1, mpc.NC, nocon);
        
        for j = 1:nuU
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('uMax%d', j);
            end
        end
        nT3Umax = length(remain);
    end
    
    conUL = [];
    nT3Umin = 0;
    if isfield(mpc, 'Umin')
        if length(mpc.Umin) ~= nuU, error('Umin dimension foul'); end
        
        [conUL remain] = addbounds(mpc.Umin, 0, mpc.NC, nocon);
        
        for j = 1:nuU
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('uMin%d', j);
            end
        end
        nT3Umin = length(remain);
    end
end
conYM = [];
nT3Ymax = 0;
conYL = [];
nT3Ymin = 0;
if nuY
    if isfield(mpc, 'Ynocon'), nocon = mpc.Ynocon;
    else nocon = [];
    end
    
    if isfield(mpc, 'Ymax') && isfield(mpc, 'Ymin')
        aux = find(mpc.Ymax < mpc.Ymin);
        if length(aux) && isfinite(mpc.Ymax(aux)) && isfinite(mpc.Ymin(aux))
            error('The bounds on Y are infeasible')
        end
    end
    if isfield(mpc, 'Ymax')
        if length(mpc.Ymax) ~= nuY, error('Ymax dimension foul'); end
        
        [conYM remain] = addbounds(mpc.Ymax, 1, mpc.OH-1, nocon); % one less Y
        
        for j = 1:nuY
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('yMax%d', j);
            end
        end
        nT3Ymax = length(remain);
    end
    
    if isfield(mpc, 'Ymin')
        if length(mpc.Ymin) ~= nuY, error('Ymin dimension foul'); end
        
        [conYL remain] = addbounds(mpc.Ymin, 0, mpc.OH-1, nocon);
        
        for j = 1:nuY
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('yMin%d', j);
            end
        end
        nT3Ymin = length(remain);
    end
end

conDUM = [];
nT3DUmax = 0;
conDUL = [];
nT3DUmin = 0;
if nudU
    if isfield(mpc, 'DUnocon'), nocon = mpc.DUnocon;
    else nocon = [];
    end
    
    if isfield(mpc, 'DUmax') && isfield(mpc, 'DUmin')
        aux = find(mpc.DUmax < mpc.DUmin);
        if length(aux) && isfinite(mpc.DUmax(aux)) && isfinite(mpc.DUmin(aux))
            error('The bounds on DU are infeasible')
        end
    end
    if isfield(mpc, 'DUmax')
        if length(mpc.DUmax) ~= nudU, error('DUmax dimension foul'); end
        
        % note: constraints on DU go only up to the control horizon (no OH)
        [conDUM remain] = addbounds(mpc.DUmax, 1, mpc.NC, nocon);
        
        for j = 1:nudU
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('DuMax%d', j);
            end
        end
        nT3DUmax = length(remain);
    end
    
    if isfield(mpc, 'DUmin')
        if length(mpc.DUmin) ~= nudU, error('DUmin dimension foul'); end
        
        [conDUL remain] = addbounds(mpc.DUmin, 0, mpc.NC, nocon);
        
        for j = 1:nudU
            if any(j == remain)
                namesT3{length(namesT3)+1} = sprintf('DuMin%d', j);
            end
        end
        nT3DUmin = length(remain);
    end
end

nuT3 = nT3Xmax+nT3Xmin+nT3Umax+nT3Umin+nT3Ymax+nT3Ymin+nT3DUmax+nT3DUmin;

% in case we have output constraint softening, expand the variable array by one (rightmost)
nslack = 0;
if length(conYL) + length(conYM) && isfield(mpc, 'YsoftP')
    if mpc.YsoftP <= 2
        % a more scalable test would be to compare it against the quadratic norm, but we don't know the objective function yet!
        error('the penalty factor YsoftP looks too small to be effective');
    elseif nT3Ymax+nT3Ymin
        disp('WARNING: did you really mean to have output parametric bounds?'); % doesn't make much sense when slackened, but what do i care?
    end
    
    nslack = 1;
    HSLACK = 1; % auxiliary quadratic term
    tmp = length(namesVars);
    namesVars{tmp+1} = 'slackY';
    % the slack doesn't affect equality constraints but we need a 0 column for it
    eqs = [eqs(:,1:tmp) zeros(ndep,nslack) eqs(:,tmp+1:end)];
end

% build all-in constraint array [variables <= RHS + T3]
% LHS columns are alinged as in the master variable index (namesVars)
constr = [];
cRHS = [];
col = 0; % count T3 thetas added
aux = (nudU+nuY)*mpc.OH - nuY; % offset to first X1
len = size(conXM, 1);
if len
    constr = [constr; zeros(len, aux) conXM(:, 1:nuX*mpc.OH) zeros(len, mpc.OH*nuU+nslack)];
    cRHS = [cRHS; [conXM(:,end) zeros(len, col) -conXM(:, nuX*mpc.OH+1:nuX*mpc.OH+nT3Xmax) zeros(len, nuT3-col-nT3Xmax)]];
    col = col + nT3Xmax;
end
len = size(conXL, 1);
if len
    constr = [constr; zeros(len, aux) conXL(:, 1:nuX*mpc.OH) zeros(len, mpc.OH*nuU+nslack)];
    cRHS = [cRHS; [conXL(:,end) zeros(len, col) -conXL(:, nuX*mpc.OH+1:nuX*mpc.OH+nT3Xmin) zeros(len, nuT3-col-nT3Xmin)]];
    col = col + nT3Xmin;
end

aux = (nudU+nuY+nuX)*mpc.OH - nuY; % offset to first U0
len = size(conUM, 1);
if len
    constr = [constr; zeros(len, aux) conUM(:, 1:nuU*mpc.NC) zeros(len, (mpc.OH-mpc.NC)*nuU+nslack)];
    cRHS = [cRHS; [conUM(:,end) zeros(len, col) -conUM(:, nuU*mpc.NC+1:nuU*mpc.NC+nT3Umax) zeros(len, nuT3-col-nT3Umax)]];
    col = col + nT3Umax;
end
len = size(conUL, 1);
if len
    constr = [constr; zeros(len, aux) conUL(:, 1:nuU*mpc.NC) zeros(len, (mpc.OH-mpc.NC)*nuU+nslack)];
    cRHS = [cRHS; [conUL(:,end) zeros(len, col) -conUL(:, nuU*mpc.NC+1:nuU*mpc.NC+nT3Umin) zeros(len, nuT3-col-nT3Umin)]];
    col = col + nT3Umin;
end

aux = nudU*mpc.OH; % offset to first Y1
len = size(conYM, 1);
if len
    % slackening output constraints: Y<=b >>> Y <= b + slack, so we need -1 coefficients for the LHS
    if nslack, tmp = -1*ones(len, nslack);
    else tmp = [];
    end
    constr = [constr; zeros(len, aux) conYM(:, 1:nuY*(mpc.OH-1)) zeros(len, mpc.OH*(nuU+nuX)) tmp];
    cRHS = [cRHS; [conYM(:,end) zeros(len, col) -conYM(:, nuY*(mpc.OH-1)+1:nuY*(mpc.OH-1)+nT3Ymax) zeros(len, nuT3-col-nT3Ymax)]];
    col = col + nT3Ymax;
end
len = size(conYL, 1);
if len
    if nslack, tmp = -1*ones(len, nslack);
    else tmp = [];
    end
    constr = [constr; zeros(len, aux) conYL(:, 1:nuY*(mpc.OH-1)) zeros(len, mpc.OH*(nuU+nuX)) tmp];
    cRHS = [cRHS; [conYL(:,end) zeros(len, col) -conYL(:, nuY*(mpc.OH-1)+1:nuY*(mpc.OH-1)+nT3Ymin) zeros(len, nuT3-col-nT3Ymin)]];
    col = col + nT3Ymin;
end

len = size(conDUM, 1);
if len
    constr = [constr; conDUM(:, 1:nuU*mpc.NC) zeros(len, (mpc.OH-mpc.NC)*nudU + (nuY+nuX+nuU)*mpc.OH-nuY+nslack)];
    cRHS = [cRHS; [conDUM(:,end) zeros(len, col) -conDUM(:, nuU*mpc.NC+1:nuU*mpc.NC+nT3DUmax) zeros(len, nuT3-col-nT3DUmax)]];
    col = col + nT3DUmax;
end
len = size(conDUL, 1);
if len
    constr = [constr; conDUL(:, 1:nuU*mpc.NC) zeros(len, (mpc.OH-mpc.NC)*nudU + (nuY+nuX+nuU)*mpc.OH-nuY+nslack)];
    cRHS = [cRHS; [conDUL(:,end) zeros(len, col) -conDUL(:, nuU*mpc.NC+1:nuU*mpc.NC+nT3DUmin) zeros(len, nuT3-col-nT3DUmin)]];
    col = col + nT3DUmin;
end

if nslack
    % slack must be positive
    % the beauty of the formulation is that we just need a single slack for all output constraints
    % it is claimed that if the penalty is big enough, when there's no constraint violation the optimal U is the same @@@ but does it work?
    constr = [constr; zeros(nslack, length(namesVars)-1) -1];
    cRHS = [cRHS; zeros(nslack, nuT3+1)];
end

if isempty(constr)
    error('Unconstrained MPC detected. Please modify');
end

%% ----------- OBJECTIVE FUNCTION -----------

namesT2 = [];
nT2Yref = 0;
nuT2 = 0;
QRf = [];
% note: mpc.QR on its own doesn't enter any setpoints
if isfield(mpc, 'Ysp_Idx') && ~isempty(mpc.Ysp_Idx)
    % tracking part: 1/2 SUM (y-ySP)*QR*(y-ySP) for i=1...OH-1
    if any(mpc.Ysp_Idx < 1) || any(mpc.Ysp_Idx > nuY) || ... % implicitly covers for nuY=0
            length(mpc.Ysp_Idx) > length(unique(mpc.Ysp_Idx))
        error('track index Ysp_Idx error');
    elseif ~isfield(mpc, 'QR') || norm(mpc.QR, 'inf') < eps
        error('SP tracking requires mpc.QR matrix');
    end
    
    % acceptable: vector of weights (for diagonal), single period matrix or all time (OH-1) matrix
    % Ysp_Idx doesn't have to be in order (sorted) but QR members should be rearranged accordingly
    nyr = length(mpc.Ysp_Idx);
    [r c] = size(mpc.QR);
    if 1==r | 1==c, QRf = diag(mpc.QR);
    elseif r ~= c | any(mpc.QR ~= mpc.QR')
        error('mpc.QR should be square and symmetric'); % symmetry imposed by transformations
    else QRf = mpc.QR;
    end
    % are we covered for all the horizon?
    aux = size(QRf, 1);
    if aux == nyr
        aux = QRf;
        
        % terminal multiplier for convenience only; same effect achieved using full length mpc.QR
        if isfield(mpc, 'pYmul'), lastQ = QRf*mpc.pYmul;
        else lastQ = QRf;
        end
        
        for i=1:mpc.OH-2
            if i < mpc.OH-2
                aux = [aux zeros(i*nyr, nyr); zeros(nyr, i*nyr) QRf];
            else
                aux = [aux zeros(i*nyr, nyr); zeros(nyr, i*nyr) lastQ];
            end
        end
        QRf = aux; % super-diagonal
    elseif aux ~= nyr*(mpc.OH-1)
        error('either one or ALL periods should be in QR');
    end
    
    % but watchout for 0s in diagonal implying a dormant (ineffective) setpoint
    for i=1:nyr
        aux = 0;
        for j=1:mpc.OH-1
            tmp = (j-1)*nyr + i;
            aux = aux + sum(abs(QRf(tmp, :))) + sum(abs(QRf(:, tmp)));
        end
        if aux < eps
            disp(mpc.Ysp_Idx(i));
            error('dormant setpoint QRi=0'); % i could drop it from Ysp_Idx instead of error?
        end
    end
    % @@@ this test isn't 100% waterproof, some setpoints may be left dangling witout QRi
    
    % expand it for all Y?
    % what if QR element zero or negative? -> not convex! e = eig(Matrices.H); if min(e) < 0
    
    % yspmode optional, fit it auto mode 1
    % for separate setpoints per interval we need the array: [{1} {2}... {N}]
    if ~isfield(mpc, 'Ysp_Mode') | isempty(mpc.Ysp_Mode)
        mpc.Ysp_Mode = [];
        for j=1:nyr
            mpc.Ysp_Mode{j} = {1:mpc.OH-1};
        end
    end
    if length(mpc.Ysp_Mode) ~= nyr
        error('Ysp_Mode invalid');
    end
    
    for j=1:nyr
        aux = [];
        for i=1:length(mpc.Ysp_Mode{j})
            aux = [aux mpc.Ysp_Mode{j}{i}];
        end
        
        if length(aux) > length(unique(aux)) | sum(aux) ~= 0.5*(mpc.OH-1)*(mpc.OH)
            disp(j);
            error('Ysp_Mode must cover all output horizon range in order'); % not really required in order!
        end
    end
    
    % values need as many members as cells in Ysp_Mode
    if ~isfield(mpc, 'Ysp_Val') | isempty(mpc.Ysp_Val)
        mpc.Ysp_Val = [];
        for j=1:nyr
            %mpc.Ysp_Val{j} = nan*[1:length(mpc.Ysp_Mode{j})]; % default to free parameters
            aux = [];
            for i=1:length(mpc.Ysp_Mode{j})
                aux = [aux {nan}];
            end
            mpc.Ysp_Val{j} = aux;
        end
    end
    
    for j=1:nyr
        if length(mpc.Ysp_Mode{j}) ~= length(mpc.Ysp_Val{j})
            disp(j);
            error('Ysp_Mode does not agree with Ysp_Val');
            % we only allow NaN (parameter) or a fixed setpoint value
            %elseif any(isinf(mpc.Ysp_Val{j})), error('nonsense in data Ysp_Val');
        end
    end
    
    % generate names for ALL theta (dependents too)
    for j=1:mpc.OH-1
        for i = mpc.Ysp_Idx
            nT2Yref = nT2Yref + 1;
            namesT2{nT2Yref} = sprintf('y%d(t+%d)_SP', i, j);
        end
    end
    
    % start generating equations for Ysp (theta-only) according to above order, including RHS
    eq_ysp = [];
    for j=1:nyr
        for i=1:length(mpc.Ysp_Mode{j})
            idx = sort(mpc.Ysp_Mode{j}{i}); % time periods to replicate first t2
            lepad = zeros(1, (idx(1)-1)*nyr + j-1);
            for k=2:length(idx)
                aux = [lepad 1 zeros(1, (idx(k)-idx(1))*nyr - 1) -1];
                eq_ysp = [eq_ysp; aux zeros(1, nT2Yref - length(aux) + 1)];
            end
            
            % if fixed, extra equation on first one (t2 = fix_value)
            if isfinite(mpc.Ysp_Val{j}{i})
                aux = [lepad 1];
                eq_ysp = [eq_ysp; aux zeros(1, nT2Yref - length(aux)) mpc.Ysp_Val{j}{i}];
            end
        end
    end
    
    % for convenience, introduce error terms for all Y-Yref
    % these are organized in time as usual but their order is like the indices
    % @@@ if you add integrators (sum of errors) they will be based on these terms
    nuE = nyr;
    % the new variables are inserted right before the (free) controls
    aux = (nudU+nuY+nuX)*mpc.OH - nuY; % offset to first U0
    eqs = [eqs(:, 1:aux) zeros(ndep, nuE*(mpc.OH-1)) eqs(:, aux+1:end)];
    % there are no constraints associated to error terms so just pad the columns
    constr = [constr(:, 1:aux) zeros(size(constr,1), nuE*(mpc.OH-1)) constr(:, aux+1:end)];
    
    tmp = [];
    for i=1:aux
        tmp{i} = namesVars{i};
    end
    for j=1:mpc.OH-1
        for i = mpc.Ysp_Idx
            tmp{length(tmp)+1} = sprintf('e%d(t+%d)', i, j);
        end
    end
    for i=aux+1:length(namesVars)
        tmp{length(tmp)+1} = namesVars{i};
    end
    namesVars = tmp;
    
    % yref parameters are now belonging to the T1 group since they are moved to the equations
    % temporarily I keep the RHS Yref coefficients in a separate matrix
    spt2LHS = zeros(ndep, nT2Yref+1); % @@@ this introduces constant RHS term! (last entry)
    lepad = zeros(nuE, nudU*mpc.OH);
    ripad = zeros(nuE, nuE*(mpc.OH-2) + nuU*mpc.OH + nslack + nuT1);
    % if Ysp_Idx was full and in order, E and Y would be identity matrices, now we have to juggle to get the Y part
    aux = eye(nuY);
    stub = -aux(mpc.Ysp_Idx, :);
    for i=1:mpc.OH-1
        eqs = [eqs; lepad stub zeros(nuE, (mpc.OH-1-i)*nuY + mpc.OH*nuX + (i-1)*nuE) eye(nuE) ripad];
        spt2LHS = [spt2LHS; zeros(nuE, (i-1)*nuE), eye(nuE) zeros(nuE,(mpc.OH-i-1)*nuE+1)];
        lepad = [lepad zeros(nuE, nuY)];
        ripad(:, 1:nuE) = [];
    end
    
    % eliminate the dependent setpoint theta, if any
    if ~isempty(eq_ysp)
        [Lsp, xid_sp, rk_sp] = U_triangulate(eq_ysp, nT2Yref, 1);
        if rk_sp < size(eq_ysp, 1) | size(eq_ysp, 2) ~= nT2Yref+1 % equations have constant term on the RHS
            error('The tracking setpoints don''t look healthy');
        end
        
        aux = spt2LHS(:, xid_sp);
        spt2LHS = [aux(:, rk_sp+1:end) - aux(:,1:rk_sp)*Lsp(:, rk_sp+1:end-1) , spt2LHS(:,end) + aux(:,1:rk_sp)*Lsp(:,end)];
        
        % lose some of the dependent names
        %namesT2{xid_sp([1:rk_sp])} = []; this doesn't shrink the cell array, so prune it manually
        tmp = [];
        aux = 1;
        % don't use setdiff because it reorganizes (sorts) the results [setdiff(1:nT2Yref, xid_sp(1:rk_sp))]
        for i=xid_sp(rk_sp+1:end)
            tmp{aux} = namesT2{i};
            aux = aux + 1;
        end
        namesT2 = tmp;
        nT2Yref = nT2Yref - rk_sp; % at the limit ALL may be fixed and gone!
    end
    
    eqs = [eqs spt2LHS]; % @@@ now with final constant term
    
    % merge old T2 theta with T1, now they appear in the equations
    for i=nuT1+1:nuT1+nT2Yref
        namesT1{i} = namesT2{i-nuT1};
    end
    nuT1 = nuT1 + nT2Yref;
    namesT2 = [];
    ndep = ndep + nuE*(mpc.OH-1);
else
    % add constant term in equations
    eqs = [eqs zeros(ndep, 1)];
    if isfield(mpc, 'QR')
        disp('WARNING: for Y setpoints you must fill in Ysp_Idx');
    end
    %@@@ begos to avoid error in line 1660!!!!
    mpc.Ysp_Idx=[];
    eq_ysp=mpc.Ysp_Idx;
end

R1f = [];
if nudU & isfield(mpc, 'R1') & norm(mpc.R1, 'inf') > eps
    % move suppression part: 1/2 SUM(i=1...OH) du*R1*du
    % note we go all the way to cover for stabilization etc. Usually for i>NC du=0 (fixed control) @@@ or perhaps this is naff and we should restrict to HC?
    % NOTE: there's no separate "rho" balancing factor; whatever you want should be premultiplied in R1
    
    [r c] = size(mpc.R1);
    if 1==r | 1==c, R1f = diag(mpc.R1);
    elseif r ~= c | any(mpc.R1 ~= mpc.R1')
        error('mpc.R1 should be square and symmetric');
    else R1f = mpc.R1;
    end
    
    % are we covered for all the horizon?
    aux = size(R1f, 1);
    if aux == nudU
        aux = R1f;
        for i=1:mpc.OH-1
            aux = [aux zeros(i*nudU, nudU); zeros(nudU, i*nudU) R1f];
        end
        R1f = aux; % super-diagonal
    elseif aux ~= nudU*mpc.OH
        error('either one or ALL periods should be in R1');
    end
end

% control term, including reference points
Rf = [];
nT2Uref = 0;
nuEU = 0;
if isfield(mpc, 'R') & norm(mpc.R, 'inf') > eps
    % move suppression part: 1/2 SUM(i=1...OH) u*R*u
    % note we go all the way to cover for stabilization etc. Usually for i>NC du=0 (fixed control) @@@ or perhaps this is naff and we should restrict to HC?
    % NOTE: there's no separate "rho" balancing factor; whatever you want should be premultiplied in R
    
    [r c] = size(mpc.R);
    if 1==r | 1==c, Rf = diag(mpc.R);
    elseif r ~= c | any(mpc.R ~= mpc.R')
        error('mpc.R should be square and symmetric');
    else Rf = mpc.R;
    end
    
    % are we covered for all the horizon?
    aux = size(Rf, 1);
    if aux == nuU
        aux = Rf;
        for i=1:mpc.OH-1
            aux = [aux zeros(i*nuU, nuU); zeros(nuU, i*nuU) Rf];
        end
        Rf = aux; % super-diagonal
    elseif aux ~= nuU*mpc.OH
        error('either one or ALL periods should be in R');
    end
    
    % control setpoints, modifying the objective as Eu'R Eu = (u-uref)R(u-uref)
    % @@@ at present, we have setpoints for ALL u and they are fixed for all the horizon
    if ~isfield(mpc, 'uRef')
        mpc.uRef = zeros(nuU, 1)';
    end
    
    if length(mpc.uRef) ~= nuU
        error('mpc.uRef should have dimension as per controls');
    end
    
    idx = find(isnan(mpc.uRef));
    
    if isempty(idx)
        idx=[];
    end
    for i = idx
        nT2Uref = nT2Uref + 1;
        namesT2{nT2Uref} = sprintf('u%d_SP', i);
    end
    
    % the control "error" terms are same in number as normal control variables (across horizon too)
    nuEU = nuU;
    
    % the new variables are inserted right before the (free) controls
    aux = (nudU+nuY+nuX)*mpc.OH - nuY + nuE*(mpc.OH-1); % offset to first U0
    eqs = [eqs(:, 1:aux) zeros(ndep, nuEU*mpc.OH) eqs(:, aux+1:end)];
    % there are no constraints associated to error terms so just pad the columns
    constr = [constr(:, 1:aux) zeros(size(constr,1), nuEU*mpc.OH) constr(:, aux+1:end)];
    
    tmp = [];
    for i=1:aux
        tmp{i} = namesVars{i};
    end
    for j=1:mpc.OH
        for i = 1:nuEU
            tmp{length(tmp)+1} = sprintf('eU%d(t+%d)', i, j-1);
        end
    end
    for i=aux+1:length(namesVars)
        tmp{length(tmp)+1} = namesVars{i};
    end
    namesVars = tmp;
    
    % like yRef, add the defining equations for EU keeping a new parameter table (uRef)
    % for the moment assume that ALL controls have a theta for their setpoint
    spt2LHS = zeros(ndep, nuEU+1); % includes a RHS term
    lepad = zeros(nuEU, aux);
    ripad = zeros(nuEU, nuU*(mpc.OH-1) + nslack + nuT1 + 1);
    for i=1:mpc.OH
        eqs = [eqs; lepad eye(nuEU) zeros(nuEU, nuEU*(mpc.OH-i) + nuU*(i-1)) -eye(nuEU) ripad];
        spt2LHS = [spt2LHS; eye(nuEU) zeros(nuEU,1)];
        lepad = [lepad zeros(nuEU, nuEU)];
        ripad(:, 1:nuEU) = [];
    end
    
    idx2 = find(~isnan(mpc.uRef));
    
    % any fixed references, update the right hand side and lose the respective column
    for i=idx2
        for j=1:mpc.OH
            spt2LHS(ndep + i + (j-1)*nuEU, end) = mpc.uRef(i);
        end
    end
    spt2LHS(:, idx2) = [];
    
    % eqs has already a constant RHS term so take care not to have TWO of them!
    tmp = eqs(:,end);
    eqs(:,end) = [];
    spt2LHS(:, end) = spt2LHS(:, end) + tmp;
    eqs = [eqs spt2LHS]; % now with final constant term
    % merge old T2 theta with T1, now they appear in the equations
    for i=nuT1+1:nuT1+nT2Uref
        namesT1{i} = namesT2{i-nuT1};
    end
    nuT1 = nuT1 + nT2Uref;
    namesT2 = [];
    ndep = ndep + nuEU*mpc.OH;
end


QXf = [];
if isfield(mpc, 'Q') & norm(mpc.Q, 'inf') > eps
    % state minimization part: 1/2 xPx + 1/2 SUM(i=1...OH-1) x*Q*x
    % this term only makes sense if 0-states are preferable, ie. we have linearized around a nominal point
    
    % acceptable: vector of weights (for diagonal), single period matrix or all time (OH-1) matrix
    [r c] = size(mpc.Q);
    if 1==r | 1==c, QXf = diag(mpc.Q);
    elseif r ~= c | any(mpc.Q ~= mpc.Q')
        error('mpc.Q should be square and symmetric');
    else QXf = mpc.Q;
    end
    
    % are we covered for all the horizon?
    aux = size(QXf, 1);
    if aux == nuX
        if ~isfield(mpc, 'P') | isempty(mpc.P)
            mpc.P = QXf;
        elseif 1==size(mpc.P,1) | 1==size(mpc.P,2)
            mpc.P = diag(mpc.P);
        end
        
        aux = QXf;
        for i=1:mpc.OH-2
            aux = [aux zeros(i*nuX, nuX); zeros(nuX, i*nuX) QXf];
        end
        if 2==mpc.OH, i=1;
        else i = i+1;
        end
        aux = [aux zeros(i*nuX, nuX); zeros(nuX, i*nuX) mpc.P];
        QXf = aux; % super-diagonal
    elseif aux ~= nuX*mpc.OH
        error('either one or ALL periods should be in Q'); % including P
    end
end


% x/Du terms are straightforward, the variables are in order as in the equations
% y terms are gatzoloptera because they (may be) less and in different order than equation vector
% from eliminating equations, all X,Y,DU and some U in the objective will be replaced by a term like {L*u + M*t1}

[L, xid, rk] = U_triangulate(eqs, size(eqs,2)-nuT1-1, 1);
if rk ~= ndep | length(xid) ~= size(eqs,2)-nuT1-1
    error('state space model or other cockup');
end%@@@begos


% all the straighforward terms e.g. xQx are extended as follows:
% (L*u+M*t1)' Q (L*u+M*t1) = u'L'QLu + t1'M'QMt1 + 2*t1'M'QLu
% so I start amassing terms for u-u, t1-t1 and t1-u

% ADDENDUM: now that equations have a constant term, we're slightly different:
% x = L*u+M*t1+C => xQx = {old} + 2*C'QLu + 2*C'QMt1 + C'QC
% so we have linear terms for u, t1 and a constant

nfree = size(eqs,2) - ndep - nuT1 -1;
if nfree <=0
    error('Something went very wrong. Contact us please.');
end
uuTerm = zeros(nfree, nfree);
ttTerm = zeros(nuT1, nuT1);
t1uTerm = zeros(nuT1, nfree); % @@@ 2 multiplier IS included
u_Lin = zeros(1, nfree);
t1_Lin = zeros(1, nuT1);
const = 0;
namesFree = [];
for i=1:nfree
    % independent variables vector may be slightly reorganized due to U blocking
    namesFree{i} = namesVars{xid(rk+i)};
end

if ~isempty(R1f)
    auxL = -L(1:nudU*mpc.OH, ndep+1:ndep+nfree);
    auxM = -L(1:nudU*mpc.OH, ndep+nfree+1:end-1);
    auxC = -L(1:nudU*mpc.OH, end);
    
    uuTerm = uuTerm + auxL'*R1f*auxL;
    ttTerm = ttTerm + auxM'*R1f*auxM;
    t1uTerm = t1uTerm + 2*auxM'*R1f*auxL;
    u_Lin = u_Lin + 2*auxC'*R1f*auxL;
    t1_Lin = t1_Lin + 2*auxC'*R1f*auxM;
    const = const + auxC'*R1f*auxC;
end

if ~isempty(QXf)
    tmp = nudU*mpc.OH + nuY*(mpc.OH-1);
    auxL = -L(tmp+1:tmp+nuX*mpc.OH, ndep+1:ndep+nfree);
    auxM = -L(tmp+1:tmp+nuX*mpc.OH, ndep+nfree+1:end-1);
    auxC = -L(tmp+1:tmp+nuX*mpc.OH, end);
    
    uuTerm = uuTerm + auxL'*QXf*auxL;
    ttTerm = ttTerm + auxM'*QXf*auxM;
    t1uTerm = t1uTerm + 2*auxM'*QXf*auxL;
    u_Lin = u_Lin + 2*auxC'*QXf*auxL;
    t1_Lin = t1_Lin + 2*auxC'*QXf*auxM;
    const = const + auxC'*QXf*auxC;
end

% now the hard bit: setpoints, expanded from (y-ySP)QR(y-ySP)
% here we have extra parameters from the T2 group
% LATEST VERSION: this is easy too now we have the error variables!
if ~isempty(QRf)
    tmp = (nudU+nuX)*mpc.OH + nuY*(mpc.OH-1);
    auxL = -L(tmp+1:tmp+nuE*(mpc.OH-1), ndep+1:ndep+nfree);
    auxM = -L(tmp+1:tmp+nuE*(mpc.OH-1), ndep+nfree+1:end-1);
    auxC = -L(tmp+1:tmp+nuE*(mpc.OH-1), end);
    
    uuTerm = uuTerm + auxL'*QRf*auxL;
    ttTerm = ttTerm + auxM'*QRf*auxM;
    t1uTerm = t1uTerm + 2*auxM'*QRf*auxL;
    u_Lin = u_Lin + 2*auxC'*QRf*auxL;
    t1_Lin = t1_Lin + 2*auxC'*QRf*auxM;
    const = const + auxC'*QRf*auxC;
end

% the control quadratic terms are now easier with the introduced EU terms!
if ~isempty(Rf)
    tmp = (nudU+nuX)*mpc.OH + (nuY+nuE)*(mpc.OH-1);
    auxL = -L(tmp+1:tmp+nuEU*mpc.OH, ndep+1:ndep+nfree);
    auxM = -L(tmp+1:tmp+nuEU*mpc.OH, ndep+nfree+1:end-1);
    auxC = -L(tmp+1:tmp+nuEU*mpc.OH, end);
    
    uuTerm = uuTerm + auxL'*Rf*auxL;
    ttTerm = ttTerm + auxM'*Rf*auxM;
    t1uTerm = t1uTerm + 2*auxM'*Rf*auxL;
    u_Lin = u_Lin + 2*auxC'*Rf*auxL;
    t1_Lin = t1_Lin + 2*auxC'*Rf*auxM;
    const = const + auxC'*Rf*auxC;
end

if norm(uuTerm, 'inf') < eps
    error('degenerate objective, did you add any Q/R matrices?');
end

% time to gather all the objective pieces in one place
% we have standard quadratic terms for u-u (uuTerm) t1-t1 (ttTerm) and t1-u (t1uTerm)
% setpoint tracking introduces t2-t2 (t2_Q) t2-u (t2uTerm) t2-t1 (t2t1Term) and linear terms t2_Lin, t1_Lin, u_Lin and t2_const (NOT!)

% we need objective of the form: 1/2 u'Hu + tFu (note 1/2 only in first part) to apply Z transformation
H = uuTerm;
if nuT2
    % merge all theta in a big array [t1 t2]
    error('these are old matrices, rearrange for u setpoints!');
    F = [t1uTerm; t2uTerm];
    TT = [ttTerm 0.5*t2t1Term'; 0.5*t2t1Term t2_Q];
    TLin = [t1_Lin t2_Lin];
    'if'
    TT
    pause
else
    F = t1uTerm;
    TT = ttTerm;
    TLin = t1_Lin;
    
    if isempty(TLin), TLin = zeros(1, nuT1); end
    
    if isempty(u_Lin), u_Lin = zeros(1, nfree); end
    if isempty(const), const = 0; end
end

if nslack
    % H as it stands now has the last row & col null, so it isn't invertible
    
    % so the only sensible alternative is to add a dummy quadratic term for it
    H(end,end) = HSLACK; 	% @@@ the magnitude of this is open to discussion!
    u_Lin(end) = mpc.YsoftP;  % linear penalty term
end

e = eig(H);
if min(e) < -ZERO
    % if we have slacks, the minimum will be zero at least!
    error('QP objective is not positive definite, watch your weights'); % QP not convex!
end

% standard QPs are in the form 1/2 xQx + fx, with the leading 1/2 coefficient implicit
% which means that we have to pass if we want to minimize xQx we have to pass 2*Q to the solver
% I have demonstrated that with both QP() and QPNAG(), passing the (x-1)^2 = x2 - 2x + 1
% qpnag(1,-2, 1, 100) -> optimal = 2 (if we pass coefficients as in the equation)
% qpnag(2*1,-2, 1, 100) -> optimal = 1 >>> the true optimal! (if we multiply the quadratic term by 2)

% premultiply with the implicit 1/2
% note that H doesn't need it because it has a ghost 1/2 in front of it
% F = 0.5*F;
TT = 0.5*TT;
TLin = 0.5*TLin;
u_Lin = 0.5*u_Lin;
const = 0.5*const;

% add the t3 vector in the picture (with zero coefficients)
TT = blkdiag(TT, zeros(nuT3, nuT3));
F = [F; zeros(nuT3, nfree)];
TLin = [TLin zeros(1,nuT3)];

nparam = nuT1 + nuT2 + nuT3;
namesTheta = namesT1;
for i=1:nuT2
    namesTheta{length(namesTheta)+1} = namesT2{i};
end
for i=1:nuT3
    namesTheta{length(namesTheta)+1} = namesT3{i};
end

% substitute dependent variables to the constraints
% Ax <= b + Ft3 => [A1 A2][xD xF]' = A1xD+A2xF <= b + Ft3, xD = LxF + Mt1

reorg = constr(:, xid); % since some U columns may be reorganized, align constraints accordingly
redA = -reorg(:,1:rk)*L(:,rk+1:rk+nfree) + reorg(:,rk+1:end);
redB = cRHS(:,1) + reorg(:,1:rk)*L(:,end);
redF = [reorg(:,1:rk)*L(:, ndep+nfree+1:end-1) zeros(length(redB), nuT2) cRHS(:, 2:end)]; % no T2 terms in constraints
% now we have Au <= b + Ft, where all theta are present and the free variables u, and NO equations

% substitutions etc may have kaked the scaling, so renormalize
tmp = [redA redF];
[rA, rb, rC, rd] = normalizeConstraints(tmp, redB, [], []);
if length(rd)
    error('we have new equations here!');
end
redA = rA(:, 1:nfree);
redF = rA(:, nfree+1:end);
redB = rb;

% store details for "online" controller (no conversion to Z)
% first inequalities Au <= b + St (no equations in sight!)
online.A = redA;
online.b = redB;
online.S = redF;
% objective function: 1/2 u'Hu + (t'F + ULin)*u + {t'TTt + Tlin*t + const}
online.H = H; % NOTE: this term is multiplied by 2!!
online.F = F;
online.ULin = u_Lin;
online.TT = TT;
online.TLin = TLin;
online.const = const;
% equations to reconstruct all variables: L [xD xF t1 c] = 0 => IxD + L1xF + L2t1 + c = 0
online.L = L;
online.xid = xid; % reordered columns in L
%online.ndep = ndep; inferred from number of equations, no?
% equations to reconstruct T2 (setpoints) Lsp [t2D t2F] = Lsp(end) -> ie includes constant RHS
if ~isempty(eq_ysp)
    online.Lsp = Lsp;
    online.xid_sp = xid_sp; % reordered T2 columns in Lsp
end
online.nuT1 = nuT1; % three groups of parameters reported
online.nuT2 = nuT2;
online.nuT3 = nuT3;
online.nParam = nparam;
online.namesTheta = namesTheta; % parameters
online.namesVars = namesVars; % all variables
online.namesFree = namesFree; % independent variables
online.mpc = mpc;
online.model = ss;

mpv = struct('c', u_Lin', 'ct', TLin', 'cc', const, 'Q', H, 'Ht', F', ...
    'Qt', TT, 'A', redA, 'b', redB, 'F', redF);

% Make the Q matrix symmetric
Q_temp = triu(mpv.Q) + tril(mpv.Q,-1)';
Q_tem = triu(Q_temp,1)/2;
mpv.Q = diag(diag(mpv.Q)) + Q_tem + Q_tem';

% For the binary case
if isfield(mpc,'BinIdx');
    Idx = [];
    for i = 1:mpc.NC
        Idx = [Idx, (i-1)*length(mpc.Umin)+mpc.BinIdx];
    end
    mpv = BinContConversion(mpv, Idx);
    
    fields = {'processing','Aeq','beq','Feq','Tmin','Tmax','original','CRA','CRb'};
    for k = 1:length(fields)
        if isfield(mpv,fields{k})
            mpv = rmfield(mpv,fields{k});
        end
    end
end
end