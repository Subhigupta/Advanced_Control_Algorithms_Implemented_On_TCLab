function [problem]=ss2qp_yalmip(mpc,syst,options)
%==========================================================================
% File Name     : <ss2qp_yalmip.m>
% Usage         : [problem]=ss2qp_yalmip(mpc,system,options)
% Description   : This function accepts an input of the MPC problem
% formulation.  It then reformulates into the mpMIQP.
%
%--------------------------------------------------------------------------
% Author        : Baris Burnak and Justin Katz
%                 Efstratios N. Pistikopoulos
% Office        : Engineering Research Building, Texas A&M University, USA
% Mail          : paroc@tamu.edu
%--------------------------------------------------------------------------
% Last Revision | Author  | Description
%---------------+---------+------------------------------------------------
% 2-Feb-2017    | JK, BB  | Initial version
% 13-Mar-2017   | JK, BB  | Minor update
% 1-Jun-2017    | JK, BB  | Output changed to problem
%==========================================================================
if nargin<3
    options=sdpsettings('solver','pop','savesolverinput',1,'verbose',1,'showprogress',1);
    options.pureexport = 1;
end

OH = mpc.OH;
NC = mpc.NC;
objective = 0;
ysp=[];
ymeas = [];
con=[];
parameters=[];
namesThita = [];
if isfield(syst,'A')
    A=syst.A;
    x = sdpvar(repmat(size(A,2),1,OH+1),repmat(1,1,OH+1));%Only really need x{1}
    Xmin = mpc.Xmin;
    Xmax = mpc.Xmax;
    if length(Xmin) ~= size(A,2)
        error('length(Xmin) does not match the number of states')
    end
    if length(Xmax) ~= size(A,2)
        error('length(Xmax) does not match the number of states')
    end
    con=[con,Xmin <= x{1} <= Xmax];
    parameters=[parameters;x{1}];
    for jj = 1:size(A,2)
        x_name = sprintf('| x%d(t) |',jj);
        namesThita = [namesThita x_name];
    end
else
    A=0;
    warning('No A matrix specified')
end

if isfield(syst,'B')
    B=syst.B;
    u = sdpvar(repmat(size(B,2),1,NC),repmat(1,1,NC));
    if NC == 1
        u_aux = cell(1,NC);
        u_aux{:} = u;
        u = u_aux;
    end
    if isfield(mpc,'Umax')
        if any(isnan(mpc.Umax))
           Umax = sdpvar(repmat(size(B,2),1,1),1); %time invariant
           nan_idx = find(isnan(mpc.Umax));
           num_idx = find(~isnan(mpc.Umax));
           Umax(num_idx) = mpc.Umax(num_idx);
           parameters = [parameters;Umax(nan_idx)];
           for jj=1:length(nan_idx)
               Umax_name = sprintf('| Umax%d |',nan_idx(jj));
               namesThita = [namesThita Umax_name];
           end
       else
           Umax = mpc.Umax;
       end
    else
        Umax = Inf*ones(size(B,2),1);
    end
    if isfield(mpc,'Umin')
        if any(isnan(mpc.Umin))
           Umin = sdpvar(repmat(size(B,2),1,1),1); %time invariant
           nan_idx = find(isnan(mpc.Umin));
           num_idx = find(~isnan(mpc.Umin));
           Umin(num_idx) = mpc.Umin(num_idx);
           parameters = [parameters;Umin(nan_idx)];
           for jj=1:length(nan_idx)
               Umin_name = sprintf('| Umin%d |',nan_idx(jj));
               namesThita = [namesThita Umin_name];
           end
       else
           Umin = mpc.Umin;
       end
    else
        Umin = -Inf*ones(size(B,2),1);
    end
    for jj = 1:NC
        con=[con,Umin <= u{jj} <= Umax];
    end
else
    B = zeros(size(A,2),1);
    warning('No B matrix specified')
end

if isfield(syst,'C')
    C=syst.C;
    if ~isfield(mpc,'DH')
        DH = 1;
        warning('Disturbance horizon is set to 1!');
    else
        DH = mpc.DH;
    end
    if OH < DH
        DH = OH;
        warning('Disturbance horizon is set to output horizon!')
    end
    d = sdpvar(repmat(size(C,2),1,DH),repmat(1,1,DH));
    if isfield(mpc,'Dmin')
        Dmin = mpc.Dmin;
    else
        Dmin = -inf*ones(size(C,2),1);
    end
    
    if isfield(mpc,'Dmax')
        Dmax = mpc.Dmax;
    else
        Dmax = inf*ones(size(C,2),1);
    end
    if DH == 1
        d_aux = cell(1,DH);
        d_aux{:} = d;
        d = d_aux;
    end
    for jj = 1:DH
        con=[con, Dmin <= d{jj} <= Dmax];
        parameters=[parameters;d{jj}];
        for kk = 1:size(C,2)
            d_name = sprintf('| d%d(t+%d) |',kk,jj-1);
            namesThita = [namesThita d_name];
        end
    end
else
    d{1}=0;
    C=zeros(size(A,2),1);
    DH = 1; %FLAG
end

if isfield(syst,'constant')
    constant = syst.constant;
else
    constant = 0;
end


if isfield(syst,'D')
    D = syst.D;
    if isfield(syst,'E')
        E = syst.E;
    else
        E = zeros(size(D,1),size(B,2));
    end
    if isfield(syst,'F')
        F = syst.F;
    else
        F = zeros(size(D,1),size(C,2));
    end
    if ~isfield(mpc,'Ymin')
        Ymin=-inf*ones(size(D,1),1);
        warning('Ymin is set to -Inf')
    else
        Ymin=mpc.Ymin;
    end
    if ~isfield(mpc,'Ymax')
        Ymax=inf*ones(size(D,1),1);
        warning('Ymax is set to Inf')
    else
        Ymax=mpc.Ymax;
    end
    if ~isfield(mpc,'QR') || ~isfield(mpc,'QRL')
        for k=1:OH-1
            dummy = 0;
            for L = 1:k
                dummy=dummy+A^(k-L)*(B*u{min(L,NC)}+C*d{min(L,DH)}+constant);
            end
            con=[con, Ymin <= D*(A^k*x{1}+dummy)+E*u{min(k+1,NC)}+F*d{min(k+1,DH)} <= Ymax ];
        end
    end
end

if isfield(mpc,'P')
    bound=OH-1;
    P=mpc.P;
    dummy=0;
    for L = 1:OH
        dummy=dummy+A^(OH-L)*(B*u{min(L,NC)}+C*d{min(L,DH)}+constant);
    end
    objective = objective + transpose(A^OH*x{1}+dummy)*P*(A^OH*x{1}+dummy);
    con=[con, Xmin<= A^OH*x{1}+dummy <= Xmax ];
else
    bound=OH;
end
for k=1:bound
    dummy = 0;
    for L = 1:k
        dummy=dummy+A^(k-L)*(B*u{min(L,NC)}+C*d{min(L,DH)}+constant);
    end
    if isfield(mpc,'Q')
        Q=mpc.Q;
        objective = objective + transpose(A^k*x{1}+dummy)*Q*(A^k*x{1}+dummy);
    end
    con=[con, Xmin<= A^k*x{1}+dummy <= Xmax ];
end

if isfield(mpc,'QR')
    if isfield(mpc,'Ymismatch')
        Yindex = 1:size(D,1);
        Ymismatch = mpc.Ymismatch;
        Ynomismatch = Yindex(~ismember(Yindex,Ymismatch));
        ymeas = sdpvar(length(Ymismatch),1);
        parameters = [parameters;ymeas];
        for jj = 1:length(Ymismatch)
            ymeas_name = sprintf('| y%d(t)_real |',Ymismatch(jj));
            namesThita = [namesThita ymeas_name];        
        end
        con = [con, Ymin(Ymismatch) <= ymeas <= Ymax(Ymismatch)];
        Yerror_mismatch = ymeas - (D(Ymismatch,:)*x{1}+E(Ymismatch,:)*u{1}+F(Ymismatch,:)*d{1});
        Yerror_nomismatch = zeros(length(Ynomismatch),1);
        Yerror = sdpvar(size(D,1),1);
        Yerror(Ymismatch) = Yerror_mismatch;
        Yerror(Ynomismatch) = Yerror_nomismatch;
    else
        Yerror = zeros(size(D,1),1);
    end
    QR=mpc.QR;
    if isfield(mpc,'Ysp_Idx')
        Ysp_Idx = mpc.Ysp_Idx;
    else
        Ysp_Idx = 1:size(D,1);
        warning('No output is specified to be tracked (mpc.Ysp_Idx=[]). All outputs are being tracked since QR matrix is provided.')
    end
    ysp=sdpvar(length(Ysp_Idx),1);
    con=[con,Ymin(Ysp_Idx) <= ysp <=Ymax(Ysp_Idx)];
    % beware: incoming patch - ON
    if isfield(syst,'E') 
        con = [con, Ymin(Ysp_Idx) <= D(Ysp_Idx,:)*x{1}+E(Ysp_Idx,:)*u{1}+F(Ysp_Idx,:)*d{min(1,DH)} <= Ymax(Ysp_Idx)];
    end
    % beware: incoming patch - OFF
    for k = 1:OH-1
        dummy=0;
        for L=1:k
            dummy=dummy+A^(k-L)*(B*u{min(L,NC)}+C*d{min(L,DH)}+constant);
        end
        ypred = D(Ysp_Idx,:)*(A^k*x{1}+dummy)+E(Ysp_Idx,:)*u{min(k+1,NC)}+F(Ysp_Idx,:)*d{min(k+1,DH)};
        yhat = ypred+Yerror(Ysp_Idx);
        objective = objective + transpose(yhat-ysp)*QR*(yhat-ysp);
        con = [con, Ymin(Ysp_Idx) <= yhat <= Ymax(Ysp_Idx)];
    end
    parameters=[parameters;ysp];
    for jj = 1:length(Ysp_Idx)
        ysp_name = sprintf('| y%d(t)_SP |',Ysp_Idx(jj));
        namesThita = [namesThita ysp_name];
    end
end

if isfield(mpc,'R')
    % not in progress anymore. needs testing though.
    if isfield(mpc,'uRef')
        if any(isnan(mpc.uRef))
            uRef = sdpvar(repmat(size(B,2),1,1),1);
        else
            uRef = mpc.uRef;
        end
        nan_idx = find(isnan(mpc.uRef));
        num_idx = find(~isnan(mpc.uRef));
        uRef(num_idx) = mpc.uRef(num_idx);
        parameters = [parameters;uRef(nan_idx)];
        con = [con, Umin(nan_idx) <= uRef(nan_idx) <= Umax(nan_idx)];
        for jj = 1:length(nan_idx)
            uRef_name = sprintf('| uRef%d(t+0) |',jj);
            namesThita = [namesThita uRef_name];
        end
    else
        uRef = 0;
    end
    R=mpc.R;
    for jj = 1:OH
        objective=objective + transpose(u{min(NC,jj)}-uRef)*R*(u{min(NC,jj)}-uRef); %Bug in ss2qp, this accounts for that by multiplying by OH
    end
end

if isfield(mpc,'R1')
% % % % % %     needs error checking. No errors encountered so far.
    R1 = mpc.R1;
    uOld = sdpvar(repmat(size(B,2),1,1),1);
    parameters = [parameters;uOld];
    con = [con, Umin <= uOld <= Umax];
    for jj = 1:size(B,2)
        uOld_name = sprintf('| uOld%d(t-1) |',jj);
        namesThita = [namesThita uOld_name];
    end
    if isfield(mpc,'DUmin')
        DUmin = mpc.DUmin;
        con = [con, DUmin <= u{1}-uOld ];
    end
    if isfield(mpc,'DUmax')
        DUmax = mpc.DUmax;
        con = [con, u{1}-uOld <= DUmax ];
    end
    objective = objective + transpose(u{1}-uOld)*R1*(u{1}-uOld);
    for jj = 2:NC
        objective = objective + transpose(u{min(NC,jj)} - u{min(NC,jj)-1})*R1*(u{min(NC,jj)} - u{min(NC,jj)-1});
        if isfield(mpc,'DUmin')
            con = [con, DUmin <= u{min(NC,jj)} - u{min(NC,jj)-1} ];
        end
        if isfield(mpc,'DUmax')
            con = [con, u{min(NC,jj)} - u{min(NC,jj)-1} <= DUmax ];
        end
    end
end

%% Linear stuff
if isfield(mpc,'QL')
    QL = mpc.QL;
    for jj = 1:OH-1
        x{jj+1} = A*x{jj} + B*u{min(jj,NC)} + C*d{min(jj,DH)} + constant;
        objective = objective + QL*x{jj+1};
        con=[con, Xmin<= x{jj+1} <= Xmax ];
    end
end

if isfield(mpc,'QRL')
    QRL = mpc.QRL;
    if isfield(mpc,'Ymismatch')
        ymeas = sdpvar(size(syst.D,1),1);
        parameters = [parameters;ymeas];
    for jj = 1:size(syst.D,1)
            ymeas_name = sprintf('| y%d(t+0)_real |',jj);
            namesThita = [namesThita ymeas_name];        
    end
        con = [con, Ymin <= ymeas <= Ymax];
        Yerror = ymeas - (D*x{1}+E*u{1}+F*d{1});
    else
        Yerror = 0;
    end
    for k = 1:OH-1
        x{k+1} = A*x{k} + B*u{min(k,NC)} + C*d{min(k,DH)} + constant;
        ypred = D*x{k+1}+E*u{min(k+1,NC)}+F*d{min(k+1,DH)};
        yhat = ypred+Yerror;
        objective = objective + QRL*(yhat);
        con = [con, Ymin <= yhat <= Ymax];
    end
end

if isfield(mpc,'RL')
    RL = mpc.RL;
    if isfield(mpc,'uRef')
        error('!uRef under construction!');
    else
        uRef = 0;
    end
    for jj = 1:NC
        objective = objective + RL*(u{min(jj,NC)} - uRef);
    end
end

if isfield(mpc,'R1L')
% % % % % %     needs error checking. No errors encountered so far.
    R1L = mpc.R1L;
    uOld = sdpvar(repmat(size(B,2),1,1),1);
    parameters = [parameters;uOld];
    for jj = 1:size(B,2)
        uOld_name = sprintf('| uOld%d(t+0) |',jj);
        namesThita = [namesThita uOld_name];
    end
    if isfield(mpc,'DUmin')
        DUmin = mpc.DUmin;
        con = [con, DUmin <= u{1}-uOld ];
    end
    if isfield(mpc,'DUmax')
        DUmax = mpc.DUmax;
        con = [con, u{1}-uOld <= DUmax ];
    end
    objective = objective + R1L*(u{1}-uOld);
    for jj = 2:OH
        objective = objective + R1L*(u{min(NC,jj)} - u{min(NC,jj)-1});
        if isfield(mpc,'DUmin')
            con = [con, DUmin <= u{min(NC,jj)} - u{min(NC,jj)-1} ];
        end
        if isfield(mpc,'DUmax')
            con = [con, u{min(NC,jj)} - u{min(NC,jj)-1} <= DUmax ];
        end
    end
end

%% Generate the optimization variables
decision = [];
for jj = 1:NC
    decision = [decision;u{jj}];
end

[~,diagnostic] = solvemp(con,objective,options,parameters,decision);

Matrices = yalmip2mpt(diagnostic);
% Some preprocessing to remove bad stuff
Matrices = removeExplorationConstraints(Matrices);
[dummy,un] = unique([Matrices.G Matrices.E Matrices.W],'rows');
Matrices.G = Matrices.G(un,:);
Matrices.E = Matrices.E(un,:);
Matrices.W = Matrices.W(un,:);
% Convert from MPT format to POP format
Matrices = mpt2pop(Matrices);

problem=Matrices;
problem.namesThita = namesThita;
end
