%==========================================================================
% File Name     : <PlotSolution.m>
% Usage         : PlotSolution(Solution, tfixed, option, online, options)
% Description   : This function plots the critical region, solution or both
% of a multi-parametric programming problem. If only a critical region is
% plotted, it does not require a solution. The options can be specified in
% the option field as follows:
%  option = {'CR','OBJ','all'}. The default is 'all'.
% The last input is also optional, and it specifies the struct 'online'
% which is obtained from the function 'ss2qp' when solving a
% multi-parametric model predictive control problem.
%
% The options are set in the function 'OptionSet.m' unless otherwise
% specified in the optional entry 'options'.
%--------------------------------------------------------------------------
% Author        : Richard Oberdieck, Nikolaos A. Diangelakis,
%                 Efstratios N. Pistikopoulos
% Office        : Engineering Research Building, Texas A&M University, USA
% Mail          : paroc@tamu.edu
%--------------------------------------------------------------------------
% Last Revision | Author  | Description
%---------------+---------+------------------------------------------------
% 02-Sep-2014   | RO      | Initial version
%---------------+---------+------------------------------------------------
% 23-Dec-2016   | NAD     | Bug fixes
%---------------+---------+------------------------------------------------
% 12-Jul-2016   | NAD     | Added feature: redundant constraints removed
%               |         | prior to CR plotting
%==========================================================================

function PlotSolution(Solution, tfixed, option, online, options)

%% Initialization
if nargin < 5
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

if nargin < 4
    online = [];
elseif isempty(option)
    option = 'all';
end

if nargin < 3
    if nargin < 2
        tfixed = NaN*ones(size(Solution(1).CR.A,2),1);
        option = 'all';
    else
        if ischar(tfixed)
            option = tfixed;
            tfixed = NaN*ones(size(Solution(1).CR.A,2),1);
        else
            option = 'all';
        end
    end
elseif isstruct(option)
    online = option;
    option = 'all';
end

% Get the dimensions of tfixed right
if size(tfixed,2) > size(tfixed,1)
    tfixed = tfixed';
end

if size(Solution(1).CR.A,2) < 2
    error('This function is only supporting problems with 2 or more parameters');
end

if size(Solution(1).CR.A,2) ~= length(tfixed)
    error('Need to specify the values of the parameters');
end

ind = find(isnan(tfixed));


%% Reduce Solution with tfixed
k = find(~isnan(tfixed)); %Finds all coefficients that are not NaN
j = isnan(tfixed); %Finds all coefficients that are NaN

for i = 1:length(Solution)
    
    Solution(i).CR.b = Solution(i).CR.b - Solution(i).CR.A(:,k) * tfixed(k);
    Solution(i).CR.A = Solution(i).CR.A(:,j);
    [Solution(i).CR.A,Solution(i).CR.b]=Redundandize(Solution(i).CR.A,Solution(i).CR.b);
    if isfield(Solution(i).CR,'P') && ~isempty(k)
        Pn = [];
        Qn = [];
        rn = [];
        for p = 1:length(Solution(i).CR.r)
            Ptemp = permute(Solution(i).CR.P(p,:,:),[2 3 1]);
            Pn = [Pn; permute(Ptemp(j,j),[3 1 2])];
            
            Qtemp = Solution(i).CR.Q(p,:);
            Qn = [Qn; Qtemp(j) + (Ptemp(k,j)+Ptemp(j,k)')*(tfixed(k)')];
            
            rtemp = Solution(i).CR.r(p);
            rn = [rn; rtemp + Qtemp(k)*tfixed(k) + tfixed(k)'*Ptemp(k,k)*tfixed(k)];
        end

        Solution(i).CR.P = Pn;
        Solution(i).CR.Q = Qn;
        Solution(i).CR.r = rn;
    end
    
    for t = 1:length(Solution(i).Solution)
        % Change X
        Solution(i).Solution(t).X(:,end) = Solution(i).Solution(t).X(:,end) + ...
            Solution(i).Solution(t).X(:,k) * tfixed(k);
        Solution(i).Solution(t).X(:,k) = [];
        
        % Change Z
        Solution(i).Solution(t).Z.d = Solution(i).Solution(t).Z.d + ...
            tfixed(k)' * Solution(i).Solution(t).Z.Q(k,k) * tfixed(k) + ...
            Solution(i).Solution(t).Z.c(k)' * tfixed(k);
        Solution(i).Solution(t).Z.c = Solution(i).Solution(t).Z.c + ...
            (Solution(i).Solution(t).Z.Q(k,:)' + ...
            Solution(i).Solution(t).Z.Q(:,k)) * tfixed(k);
        Solution(i).Solution(t).Z.c(k) = [];
        Solution(i).Solution(t).Z.Q = Solution(i).Solution(t).Z.Q(j,j);
    end
end

%% Define function Z
Z = @(Q,c,d,t) t' * Q * t + c' * t + d;

%% Treat Infeasibilities
% Infeasible regions are stored in InfesCR
q = 1;
InfesCR = [];
while q <= length(Solution)
    % If it is infeasible (inf-values), then we only have one solution
    
    if any(any(isinf(Solution(q).Solution(1).X)))
        if isempty(InfesCR);
            InfesCR = Solution(q);
        else
            InfesCR(end+1) = Solution(q);
        end
        Solution(q) = [];
        q = q - 1;
    end
    q = q + 1;
end
idx_max = num2str(numel(num2str(length(Solution))));

%% Colours and Start (copied from Martina)
% instead of using just the 'ymcrgbkw' colors, go through the HSV instead
% one color for each CR
colors = [];
aux = hsv;
for i=1:length(Solution)
    % pseudo-random color assignment idea direct from switzerlanden :)
    j = mod(i*7, size(aux,1)) + 1;
    colors(i,:) = aux(j,:);
end

% disp('Creating solution figure...');
figure; % make sure a fresh one is created

legLab = [];
q = 1;

%% For each critical region
if sum(isnan(tfixed)) == 1
    q = 1;
    handl = [];
    for i = 1:length(Solution)
        
        [mi,ma] = Const2Bounds(Solution(i).CR.A,Solution(i).CR.b,options);
        cra = zeros(100,1);
        if strcmp(option,'OBJ') || strcmp(option,'all') || isfield(Solution(i).CR,'P')
            pts = linspace(mi,ma,100);
            zval = zeros(100,1);
            xval = zeros(100,4);
            for u = 1:length(pts)
                [~,x1,z1] = PointLocation(Solution(i),pts(u),options);
                if isempty(z1)
                    cra(u) = NaN;
                    zval(u) = NaN;
                    xval(u,:) = NaN;
                else
                    zval(u) = z1;
                    xval(u,:) = x1';
                end
            end
        end
        
        if ~any(isinf([mi;ma])) && any(~isnan(cra))
            if strcmp(option,'CR') || strcmp(option,'all')
                % Only create a subplot if option is set to 'all'
                if strcmp(option,'all')
                    subplot(1,2,1)
                end
                
                if isfield(Solution(i).CR,'P')
                    ii = find(~isnan(cra),1,'first');
                    ij = find(~isnan(cra),1,'last');
                    plot(pts(ii),cra(ii),'o-','Color',colors(i,:),'LineWidth',2,'MarkerFaceColor',colors(i,:));
                    hold on
                    plot(pts(ij),cra(ij),'o-','Color',colors(i,:),'LineWidth',2,'MarkerFaceColor',colors(i,:));
                    plot(pts,cra,'Color',colors(i,:),'LineWidth',2,'MarkerFaceColor',colors(i,:));
                else
                    plot([mi;ma],[0;0],'o-','Color',colors(i,:),'LineWidth',2,'MarkerFaceColor',colors(i,:));
                end
                hold on
                q = q + 1;
            end
            
            
            if strcmp(option,'OBJ') || strcmp(option,'all')
                if strcmp(option,'all')
                    subplot(1,2,2)
                end

                %plot(pts,zval,'Color',colors(i,:),'LineWidth',2);
                plot(pts,xval(:,2),'b','LineWidth',2);
                keyboard
                hold on
            end
            legLab = [legLab; sprintf(['CR%0',idx_max,'d'], i)]; % Name the critical region
        end
    end
    
    legend(legLab);
    if strcmp(option,'CR') || strcmp(option,'all')
        if strcmp(option,'all')
            subplot(1,2,1)
        end
        
        if isempty(online)
            t1 = num2str(ind(1));
            t1 = ['\theta_{',t1,'}'];
            xlabel(t1,'Interpreter','tex');
        else
            t1 = online.namesTheta{ind(1)};
            xlabel(t1,'Interpreter','None');
        end
        ylabel(' ');
        set(gca,'ytick',[]);
        title(sprintf('%d Feasible Region Fragments', q-1));
        grid on;
    end
    
    if strcmp(option,'OBJ') || strcmp(option,'all')
        if strcmp(option,'all')
            subplot(1,2,2)
        end
        if isempty(online)
            t1 = num2str(ind(1));
            t1 = ['\theta_{',t1,'}'];
            xlabel(t1,'Interpreter','tex');
        else
            t1 = online.namesTheta{ind(1)};
            xlabel(t1,'Interpreter','None');
        end
        ylabel('z(\theta)');
        title(sprintf('Objective Function Value'));
    end
    
    
    if any(~isnan(tfixed))
        str = 'Fixed:';
        k = find(~isnan(tfixed));
        for j = 1:length(k)
            s1 = num2str(k(j));
            s2 = num2str(tfixed(k(j)));
            if isempty(online)
                str = [str; {['\theta_{',s1,'}: {',s2,'}']}];
            else
                str = [str; {[online.namesTheta{k(j)},': ',s2]}];
            end
        end
        if isempty(online)
            annotation('textbox', [0.2,0.4,0.1,0.1],...
                'String', str, 'BackgroundColor','white','Interpreter','tex');
        else
            annotation('textbox', [0.2,0.4,0.1,0.1],...
                'String', str, 'BackgroundColor','white','Interpreter','none');
        end
        plotedit('on');
    end
    hold off
elseif sum(isnan(tfixed)) == 2
    for i = 1:length(Solution)
        CR.A = Solution(i).CR.A;
        CR.b = Solution(i).CR.b;
        exact = false;
        if isfield(Solution(i).CR,'P')
            exact = true;
        end
        
        % Check whether the CR is empty
        [~,r] = Chebyshev(CR.A, CR.b,[],[],options);
        
        if ~isnan(r) && r > options.tolerance
            %% Preparation of plotting
            % Find vertices of the critical region
            try
                V = con2vert(CR.A, CR.b);
                
                % Extract the x and y values from V
                px = V(:,1)';
                py = V(:,2)';
                
                % Get the convex hull of the respective points. This does nothing else
                % than order the vertices in a 'counter-clockwise' manner.
                nop = false;
                K = convhull(px, py);
            catch
                % Sometimes con2vert dies, which is why we have the 'try-catch'
                % routine
                nop = true;
            end
            if ~nop
                px = px(K);
                py = py(K);
                
                % Create Grid for surf
                Max_Theta = max([px; py],[],2);
                Min_Theta = min([px; py],[],2);
                
                % Create Mesh
                points = 100;
                Vector_X = linspace(Min_Theta(1), Max_Theta(1), points);
                Vector_Y = linspace(Min_Theta(2), Max_Theta(2), points);
                [Tx, Ty] = meshgrid(Vector_X, Vector_Y); %Square grid
                [aux1, aux2] = size(Tx);
                in = inpolygon(Tx,Ty,px',py');
                B = NaN * ones(points);
                
                %% Obtain optimal objective function value
                Z_matrix = inf * ones(aux1,aux2);
                if strcmp(option,'OBJ') || strcmp(option,'all') || exact
                    X = Tx;
                    Y = Ty;
                    for k = 1:aux1
                        for l = 1:aux2
                            if isfield(Solution(i).CR,'P')
                                in_it = (all(Solution(i).CR.A*[Tx(k,l);Ty(k,l)] - Solution(i).CR.b ...
                                    <= options.tolerance) && all(exactLocate([Tx(k,l);Ty(k,l)], ...
                                    Solution(q).CR.P,Solution(i).CR.Q,Solution(i).CR.r) ...
                                    <= options.tolerance));
                            else
                                in_it = all(Solution(q).CR.A*[Tx(k,l);Ty(k,l)] <= Solution(q).CR.b + options.tolerance);
                            end
                            
                            if in_it
                                if strcmp(option,'OBJ') || strcmp(option,'all')
                                    for ku = 1:length(Solution(i).Solution)
                                        Q = Solution(i).Solution(ku).Z.Q;
                                        c = Solution(i).Solution(ku).Z.c;
                                        d = Solution(i).Solution(ku).Z.d;
                                        
                                        Z_matrix(k,l) = min([Z_matrix(k,l) Z(Q, c, d, [Tx(k,l);Ty(k,l)])]);
                                    end
                                end
                            else
                                Z_matrix(k,l) = NaN;
                                X(k,l) = NaN;
                                Y(k,l) = NaN;
                                in(k,l) = false;
                            end
                        end
                    end
                    X = X(~isnan(X));
                    Y = Y(~isnan(Y));
                end
                
                %% Plot the CR (distinguish between exact and not due to affinity)
                if strcmp(option,'CR') || strcmp(option,'all')
                    % Only create a subplot if option is set to 'all'
                    if strcmp(option,'all')
                        subplot(1,2,1)
                    end
                    if exact
                        if ~all(all(in == 0))
                            B(in) = 1;
                            mesh(Tx, Ty, B, 'FaceColor', colors(i,:), 'EdgeColor','none')
                        else
                            q = q - 1;
                        end
                    else
                        fill(px, py, colors(i,:)); % instead of color name, use HSV vector equivalent
                    end
                    legLab = [legLab; sprintf(['CR%0',idx_max,'d'], i)]; % Name the critical region
                end
                hold on;
                q = q + 1;
            end
            %% Plot the objective function
            if strcmp(option,'OBJ') || strcmp(option,'all')
                if strcmp(option,'all')
                    subplot(1,2,2)
                end
                surf(Tx,Ty,Z_matrix,'EdgeColor','none');
                hold on;
            end
        end
    end
    
    %% Post-Processing of plots
    if strcmp(option,'CR') || strcmp(option,'all')
        if strcmp(option,'all')
            subplot(1,2,1)
        end
        if isfield(Solution,'CR_Exact')
            az = 0;
            el = 90;
            view(az, el);
        end
        
        legend(legLab);
        if isempty(online)
            t1 = num2str(ind(1));
            t1 = ['\theta_{',t1,'}'];
            t2 = num2str(ind(2));
            t2 = ['\theta_{',t2,'}'];
            xlabel(t1,'Interpreter','tex');
            ylabel(t2,'Interpreter','tex');
        else
            t1 = online.namesTheta{ind(1)};
            t2 = online.namesTheta{ind(2)};
            xlabel(t1,'Interpreter','None');
            ylabel(t2,'Interpreter','None');
        end
        
        
        title(sprintf('%d Feasible Region Fragments', q-1));
        grid on;
    end
    
    if strcmp(option,'OBJ') || strcmp(option,'all')
        if strcmp(option,'all')
            subplot(1,2,2)
        end
        if isempty(online)
            t1 = num2str(ind(1));
            t1 = ['\theta_{',t1,'}'];
            t2 = num2str(ind(2));
            t2 = ['\theta_{',t2,'}'];
            xlabel(t1,'Interpreter','tex');
            ylabel(t2,'Interpreter','tex');
        else
            t1 = online.namesTheta{ind(1)};
            t2 = online.namesTheta{ind(2)};
            xlabel(t1,'Interpreter','None');
            ylabel(t2,'Interpreter','None');
        end
        zlabel('z(\theta)');
        title(sprintf('Objective Function Value'));
    end
    
    if any(~isnan(tfixed))
        str = 'Fixed:';
        k = find(~isnan(tfixed));
        for j = 1:length(k)
            s1 = num2str(k(j));
            s2 = num2str(tfixed(k(j)));
            if isempty(online)
                str = [str; {['\theta_{',s1,'}: {',s2,'}']}];
            else
                str = [str; {[online.namesTheta{k(j)},': ',s2]}];
            end
        end
        if isempty(online)
            annotation('textbox', [0.2,0.4,0.1,0.1],...
                'String', str, 'BackgroundColor','white','Interpreter','tex');
        else
            annotation('textbox', [0.2,0.4,0.1,0.1],...
                'String', str, 'BackgroundColor','white','Interpreter','none');
        end
        plotedit('on');
    end
elseif sum(isnan(tfixed)) == 3
    if strcmp(option,'OBJ')
        disp('Objective function plot not possible for 3 parameters');
        return
    end
    if isfield(Solution,'CR_Exact')
        disp('The plot of the exact solution for 3 parameters is currently not supported.');
        return
    end
    q = 1;
    handl = [];
    for i = 1:length(Solution)
        
        XH.A = Solution(i).CR.A;
        XH.b = Solution(i).CR.b;
        
        V = con2vert(XH.A,XH.b);
        
        % Match indices with half-spaces
        idx = {1:2, [1,3], 2:3, 1:3};
        for k = 1:size(XH.A,1)
            
            VX = V(abs(XH.A(k,:)*V'-XH.b(k)) <= 1e-3,:);
            if ~isempty(VX)
                for j = 4:-1:1
                    try
                        [K,Vc] = convhull(VX(:,idx{j}));
                        if Vc < 1e-7
                            continue
                        end
                        break
                    catch
                        continue
                    end
                end
                if Vc > 1e-7
                    if (max(VX(K,2)) - min(VX(K,2)) > 1e-3)
                        hp = fill3(VX(K,1), VX(K,2), VX(K,3), colors(i,:));
                        alpha(hp,0.5);
                        if k == 1
                            handl(i) = hp;
                        end
                    end
                end
            end
            hold on
        end
        legLab = [legLab; sprintf(['CR%0',idx_max,'d'], i)]; % Name the critical region
        q = q + 1;
    end
    
    legend(handl,legLab);
    if isempty(online)
        t1 = num2str(ind(1));
        t1 = ['\theta_{',t1,'}'];
        t2 = num2str(ind(2));
        t2 = ['\theta_{',t2,'}'];
        t3 = num2str(ind(3));
        t3 = ['\theta_{',t3,'}'];
        xlabel(t1,'Interpreter','tex');
        ylabel(t2,'Interpreter','tex');
        zlabel(t3,'Interpreter','tex');
    else
        t1 = online.namesTheta{ind(1)};
        t2 = online.namesTheta{ind(2)};
        t3 = online.namesTheta{ind(3)};
        xlabel(t1,'Interpreter','None');
        ylabel(t2,'Interpreter','None');
        zlabel(t3,'Interpreter','None');
    end
    
    title(sprintf('%d Feasible Region Fragments', q-1));
    grid on;
    
    if any(~isnan(tfixed))
        str = 'Fixed:';
        k = find(~isnan(tfixed));
        for j = 1:length(k)
            s1 = num2str(k(j));
            s2 = num2str(tfixed(k(j)));
            if isempty(online)
                str = [str; {['\theta_{',s1,'}: {',s2,'}']}];
            else
                str = [str; {[online.namesTheta{k(j)},': ',s2]}];
            end
        end
        if isempty(online)
            annotation('textbox', [0.2,0.4,0.1,0.1],...
                'String', str, 'BackgroundColor','white','Interpreter','tex');
        else
            annotation('textbox', [0.2,0.4,0.1,0.1],...
                'String', str, 'BackgroundColor','white','Interpreter','none');
        end
        plotedit('on');
    end
else
    disp('There are too many variable parameters, please fix using ''tfixed''.');
end
end