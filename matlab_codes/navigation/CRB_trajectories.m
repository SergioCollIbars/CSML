clear;
clc;
close all;
addpath("data/")
addpath("functions/")
addpath("../../QGG_navigation/functions/measurements/")
addpath("../../QGG_navigation/functions/")
addpath("../../QGG_navigation/data/")
set(0,'defaultAxesFontSize',16);
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
%%                        PCRLB UPPER BOUND
% Description: Test the sequential CRLB in different set of trajectories
% download from JPL website.
% Author: Sergio Coll Ibars
% Date: 02/17/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization units
GM_E = 3.986004418E14;              % [m^3/s^2]
GM_M = 4.9048695E12;                % [m^3/s^2]
D    = 384399E3;                    % [m]
n    = sqrt((GM_E + GM_M) / D^3);   % [1/s]
Re = 6371E3;                        % Earth radii[m]
Re_ND = Re / D;                     % [-]
Rm = 1740E3;                        % Moon radii [m]
Rm_ND = Rm / D;                     % [-]

[planetParams, poleParams, Cmat, Smat, TIME, ~] = ...
    load_universe("CR3BP", [0, pi], 1);

% load initial condition and trajectory
data = load('EM_NHalo_south_L2_Family.mat');
% % data = load('JPL_EM_Lyap_L1_Family.mat');
% % data = load('JPL_EM_Lyap_L2_Family.mat');
% % data = load('JPL_EM_Vert_L2_Family.mat');

Nd =length(data.trajFam);
index = 1:1:Nd;
index = 1:40:Nd;
% % index = [1, 100, 500, Nd];
periapsis = 1;  % starting @ periapsis? 1 yes / 0 no
insidePlanet = zeros(1, length(index));

% trajectories struct
dataOrbit = struct('traj', nan(3, 1E4), 'time', nan(1, 1E4), 'JC', []);
dataCRB_pos    = struct('value', nan(1, 1E4));
period_mat     = ones(1, length(index));
for i = 1:length(index)
    % Filling with NaN data 
    dataOrbit(i).traj = ones(3, 1E4)*NaN;
    dataOrbit(i).time = ones(1, 1E4)*NaN;
    dataOrbit(i).JC   = data.trajFam{index(i), 1}.CJ; 
    dataCRB_pos(i).value   = ones(1, 1E4) * NaN;
    period_mat(i) = data.trajFam{index(i), 1}.period; 
end
maxPeriod = max(period_mat); minPeriod = min(period_mat);
dataCRB_vel = dataCRB_pos;

for j = 1:length(index)
    disp('Progress ... ' + string(j/length(index) * 100));

    x0   = data.trajFam{index(j),1}.iState;
    mu = data.trajFam{index(j),1}.mu;
    P    = data.trajFam{index(j), 1}.period;
    Ns   = 6;
    planetParams(1) = mu;
    
    [r0, v0] = rotate2inertial(x0(1:3), x0(4:6), 0, 1);
    X0 = [r0;v0];
    
    if(periapsis == 0)
        % integrate trajectory. Take the state @ periapsis
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        system = "CR3BP";
        STM0 = reshape(eye(6,6), [36, 1]);
        [t, state] = ode113(@(t, x) EOM_3BP(t, x, planetParams, ...
            poleParams, Cmat, Smat, system), [0, 1*P], [x0(1:6); STM0], options);
        
        % find minum position w.r.t Moon
        [~, posM] = compute_EM_position(0, mu);
        relPos_M = vecnorm(state(:, 1:3)' - posM);
        [~, idx] = min(relPos_M);

        % get the new initial conditions. Inertial coordinates
        [r0, v0] = rotate2inertial(state(idx, 1:3)', state(idx, 4:6)', t(idx), 1);
        X0 = [r0; v0];
        tmin = t(idx);   % [-]
    else
        tmin = 0;                % [-]
    end
    
    % Simulation parameters
    tmax = 2*maxPeriod + tmin;         % [-]
% %     tmax = 1.5*P;
    meas = "QGG";            % QGG / DSN
    
    % integrate trajectory.
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    system = "CR3BP_inertial";
    STM0 = reshape(eye(6,6), [36, 1]);
    [t, state] = ode113(@(t, x) EOM_3BP(t, x, planetParams, ...
        poleParams, Cmat, Smat, system), [tmin, tmax], [X0; STM0], options);
    TIME = t;
    dataOrbit(j).traj = state(:, 1:3)';
    dataOrbit(j).time = t';
    
    % compute Earth and Moon positions in inertial frame
    [posE, posM] = compute_EM_position(t, mu);

    relPos_E = vecnorm(state(:, 1:3)' - posE);
    relPos_M = vecnorm(state(:, 1:3)' - posM);

    
    if(sum(relPos_E < Re_ND) || sum(relPos_M < Rm_ND))
        insidePlanet(j) = 1;
    end
    
    if(insidePlanet(j) == 0)
        % initial uncetainty and meas weight
        sG    = 1E-12 / (n^2);                         % [-]
        R_QGG = diag([sG, sG, sG, sG, sG, sG].^2);     % [-]
        
        sP = 1E6/D;                                   % [-]
        sV = 10/(D*n);                                % [-]
        P0 = diag([sP, sP, sP, sV, sV, sV].^2);       % [-]
        
        % compute Upper bounds
        scaleP = D;     % [m]
        scaleV = D*n;   % [m/s]
        [CRB_pos, CRB_vel] = compute_CRB_Linear(state, posM, posE, t, P0, R_QGG, Ns,...
            mu, meas);
% %         [X, P, CRB_pos, CRB_vel] = compute_CRB_Unscented(X0, P0, ...
% %             planetParams, t, R_QGG);
        dataCRB_pos(j).value = CRB_pos.*scaleP;    % [m]
        dataCRB_vel(j).value = CRB_vel.*scaleV;    % [m/s]
    else
        dataOrbit(j).traj = ones(3, length(t)) * NaN;
    end
end


%%  plot trajectory

figure()
colormap("jet");
for j = 1:length(index)
    if(insidePlanet(j) == 0)
        Nt = length(dataOrbit(j).time);
        
        % get trajectory
        [rB] = rotate2BodyFrame(dataOrbit(j).time,  dataOrbit(j).traj');
        rB = rB.*(scaleP./1E3);                             % [km]
        
        % Detect last orbital period
        T = period_mat(j);
        time = dataOrbit(j).time./T;
        maxP = max(time);
        minP = maxP - 1;
        [~, idx] = min(abs(time - minP));  % index of closest element

        plot3(rB(1, idx:end), rB(2, idx:end), rB(3, idx:end), ...
            'Color', 'k');
        axis equal;
        hold on;
    end
end
plot(-mu*(scaleP./1E3), 0, "o",'MarkerFaceColor',"#7E2F8E", 'MarkerEdgeColor', "#7E2F8E")
plot((1-mu)*(scaleP./1E3),0, "o",'MarkerFaceColor',"#77AC30", 'MarkerEdgeColor', "#77AC30")
title('Orbit family')
xlabel('X [Km]')
ylabel('Y [Km]')

figure()
colormap("jet");
maxVal = ones(1, length(index)) * NaN;
minVal = maxVal;
for j = 1:length(index)
    if(insidePlanet(j) == 0)
        Nt = length(dataOrbit(j).time);
        
        % get trajectory
        [rB] = rotate2BodyFrame(dataOrbit(j).time,  dataOrbit(j).traj');
        rB = rB.*(scaleP./1E3);                             % [km]
        scalarValues = dataCRB_pos(j).value;
        maxVal(j) = max(scalarValues);
        minVal(j) = min(scalarValues);
        
        % Detect last orbital period
        T = period_mat(j);
        time = dataOrbit(j).time./T;
        maxP = max(time);
        minP = maxP - 1;
        [~, idx] = min(abs(time - minP));  % index of closest element

        scatter3(rB(1, idx:end), rB(2, idx:end), rB(3, idx:end), 7, log10(scalarValues(idx:end)), 'filled');
        axis equal;
        hold on;
    end
end
% % plot(-mu*(scaleP./1E3), 0, "o",'MarkerFaceColor',"#7E2F8E", 'MarkerEdgeColor', "#7E2F8E")
plot((1-mu)*(scaleP./1E3),0, "o",'MarkerFaceColor',"#77AC30", 'MarkerEdgeColor', "#77AC30")
h= colorbar; % Show color scale
clim([-2, 6]);
h.Ticks = linspace(min(clim), max(clim), 9); % Set tick positions
h.TickLabels = {'< 1 cm', '10 cm', '1m', '10 m', '100 m', '1 km', '10 km', '100 km', '+ 1000 km'}; % Custom labels
title('Maximum position value from covariance')
xlabel('X [Km]')
ylabel('Y [Km]')



% plot trajectory
figure()
colormap("jet");
maxVal = ones(1, length(index)) * NaN;
minVal = maxVal;
for j = 1:length(index)
    if(insidePlanet(j) == 0)
        Nt = length(dataOrbit(j).time);
        
        % get trajectory
        [rB] = rotate2BodyFrame(dataOrbit(j).time,  dataOrbit(j).traj');
        rB = rB.*(scaleP./1E3);                             % [km]
        scalarValues = dataCRB_vel(j).value;
        maxVal(j) = max(scalarValues);
        minVal(j) = min(scalarValues);

        % Detect last orbital period
        T = period_mat(j);
        time = dataOrbit(j).time./T;
        maxP = max(time);
        minP = maxP - 1;
        [~, idx] = min(abs(time - minP));  % index of closest element

        scatter3(rB(1, idx:end), rB(2, idx:end), rB(3, idx:end), ...
            7, log10(scalarValues(idx:end)), 'filled');
        axis equal;
        hold on;
    end
end
% % plot(-mu*(scaleP./1E3), 0, "o",'MarkerFaceColor',"#7E2F8E", 'MarkerEdgeColor', "#7E2F8E")
plot((1-mu)*(scaleP./1E3),0, "o",'MarkerFaceColor',"#77AC30", 'MarkerEdgeColor', "#77AC30")
h= colorbar; % Show color scale
clim([-3, 1]);
h.Ticks = linspace(min(clim), max(clim), 5); % Set tick positions
h.TickLabels = {'< 1 mm/s', '1 cm/s', '10 cm/s', '1 m /s', '10 m /s'}; % Custom labels
xlabel('X [Km]')
ylabel('Y [Km]')
title('Maximum velocity value from covariance')

% Plot position accuracy over time
figure
scalarVals = ones(1, length(index)) * NaN;
for j = 1:length(index)
    scalarVals(j) = dataOrbit(j).JC;
end
cmap = turbo;                      % or parula, turbo, etc.
nColors = size(cmap,1);

% Normalize scalar array to [1, nColors]
cmin = min(scalarVals);
cmax = max(scalarVals);
for j = 1:length(index)
    if insidePlanet(j)==0
        ax(1) = subplot(2, 1, 1);
        % Map scalar value to a row of the colormap
        cIdx = round( 1 + (scalarVals(j)-cmin) * (nColors-1) / (cmax-cmin) );
        cIdx = max(min(cIdx,nColors),1);  % safety
        thisColor = cmap(cIdx,:);
        timeX  = dataOrbit(j).time./planetParams(3) / 86400;
        semilogy(timeX, dataCRB_pos(j).value, ...
                 'LineWidth',2,'Color',thisColor);
        hold on
    end
end
ylabel('[m]')
xlabel('Time [days]')

for j = 1:length(index)
    if insidePlanet(j)==0
        ax(2) = subplot(2, 1, 2);
        % Map scalar value to a row of the colormap
        cIdx = round( 1 + (scalarVals(j)-cmin) * (nColors-1) / (cmax-cmin) );
        cIdx = max(min(cIdx,nColors),1);  % safety
        thisColor = cmap(cIdx,:);
        timeX  = dataOrbit(j).time./planetParams(3) / 86400;
        semilogy(timeX, dataCRB_vel(j).value, ...
                 'LineWidth',2,'Color',thisColor);
        hold on
    end
end
ylabel('[m/s]')
xlabel('Time [days]')
colormap(cmap)                       % set figure colormap
cb = colorbar;                       % add colorbar
cb.Position = [0.92 0.11 0.02 0.77]; % adjust numbers to span both
clim([cmin cmax])                    % match colorbar limits to data
sgtitle('Position and velocity uncertainty')

if(length(index) < 100)
    figure()
    for j = 1:length(index)
        if(insidePlanet(j) == 0)
            Nt = length(dataOrbit(j).time);
            [rB] = rotate2BodyFrame(dataOrbit(j).time,  dataOrbit(j).traj');
            rB = rB.*(scaleP./1E3);                             % [km]
            cIdx = round( 1 + (scalarVals(j)-cmin) * (nColors-1) / (cmax-cmin) );
            cIdx = max(min(cIdx,nColors),1);  % safety
            thisColor = cmap(cIdx,:);
            plot3(rB(1, :), rB(2, :), rB(3, :), 'Color', thisColor, ...
                'LineWidth', 1.5)
            axis equal;
            hold on;
        end
    end
    grid on;
    xlabel('X [km]'); ylabel('Y [Km]'); zlabel('[Km]')
    title('3D orbit');
    plot((1-mu)*(scaleP./1E3),0, "o",'MarkerFaceColor',"#77AC30",...
        'MarkerEdgeColor', "#77AC30")   % plot Moon
    plot(-mu*(scaleP./1E3), 0, "o",'MarkerFaceColor',"#7E2F8E", ...
        'MarkerEdgeColor', "#7E2F8E")   % plot Earth
end

%% FUNCTIONS
function [rB] = rotate2BodyFrame(t, state)
    Nt = length(t);
    rB = ones(3, Nt) * NaN;

    for k =1:Nt
        M = t(k);
        NB = [cos(M), -sin(M), 0;...
            sin(M), cos(M), 0;...
            0, 0, 1];
        rB(:, k) = NB' * state(k, 1:3)'; % [-]
    end
    
end

function [posE, posM] = compute_EM_position(t, mu)
    Nt = length(t);
    posE = ones(3, Nt) * NaN;
    posM = ones(3, Nt) * NaN;

    for k = 1:Nt
        M = t(k);  % [rad]

        % Earth position
        posE(:, k) = -[mu*cos(M);mu*sin(M);0];
    
        % Moon position
        posM(:, k)  = -[(mu-1)*cos(M);(mu-1)*sin(M); 0];
    end
end

function [CRB_pos, CRB_vel] = compute_CRB_Linear(state, posM, posE, t, P0, R0, Ns, mu, meas)
    % compute PCRLB and upper bound to compare
    P0 = P0(1:Ns, 1:Ns);
    Nt  = length(t);
    STM = state(:, 7:end);
    pos = state(:, 1:3)';
    vel = state(:, 4:6)';
    PCRLB = zeros(Nt, Ns*Ns); PCRLB(1, :) = reshape(inv(P0), [Ns*Ns, 1]);
    CRB_pos    = zeros(Nt, 1); CRB_vel = CRB_pos;
    for k = 1:Nt
         if(k == 1)
            PHI_1 = eye(6,6);
            Aiprev_plus = reshape(PCRLB(k, :), [Ns,Ns]);
        else
            PHI_1 = reshape(STM(k-1, :), [6,6]);            % from 0 to k-1
            Aiprev_plus = reshape(PCRLB(k-1, :), [Ns,Ns]);
         end
        PHI_2 = reshape(STM(k, :), [6, 6]);     % from 0 to k
        PHI = PHI_2 / (PHI_1);               % from k-1 to k
        PHI = PHI(1:Ns, 1:Ns);
        PHI_inv = inv(PHI);                     % from k to k-1
    
        % measurement partials
        if(meas == "QGG")
            posRel = pos(:, k) - posE(:, k); % relative position to Earth
            h1 = compute_grad_posPartials(1-mu, posRel(1), posRel(2), posRel(3));
    
            posRel = pos(:, k) - posM(:, k); % relative position to Moon
            h2 = compute_grad_posPartials(mu, posRel(1), posRel(2), posRel(3));
            if(Ns == 6), h = [h1 + h2, zeros(6, 3)]; else, h = h1+h2; end
        elseif(meas == "DSN")
            posRel = pos(:, k) - [-mu;0;0];
            velRel = vel(:, k);
            h = computePartials_DSN(posRel, velRel);
            h = h(1:2, 1:Ns); 
        end
        
        % IF filter
        Ai_min = PHI_inv' * Aiprev_plus * PHI_inv;
        Ai_plus = Ai_min + h' * (R0 \ h);
    
        p = diag(inv(Ai_plus));

        CRB_pos(k)  = max(p(1:3))^(1/2);       % max direction constraint   [m]
        if(Ns == 6)
            CRB_vel(k)  = max(p(4:6))^(1/2);   % max direction constraint   [m/s]
        else
            CRB_vel(k)  = 1;                   % max direction constraint   [m/s]
        end
    
        % update information matrix
        PCRLB(k, :)=  reshape(Ai_plus, [Ns*Ns, 1]);
    end
end

function [X, P, CRB_pos, CRB_vel] = compute_CRB_Unscented(X0, P0, planetParams, t, R0)
    % compute PCRLB and upper bound to compare
    Nt  = length(t);
    Ns  = length(X0);
    Nm  = 6;
    
    % output values
    X = zeros(6, Nt); X(:, 1) = X0;
    P = zeros(Nt, 6*6); P(1, :) = reshape(P0, [36, 1]);
    CRB_pos = ones(1, Nt) * NaN;
    CRB_vel = CRB_pos;

    for k = 2:Nt
        % create sigma points. State @ k-1
        Xhat_prev = X(:, k-1);
        [xhat_i_prev] = sigmaPoints_state(Ns, P(k-1, :), Xhat_prev);
    
        % propagate sigma points using N.L funcitons. State @ k
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        time = [t(k-1), t(k)];
        Xhat_i = xhat_i_prev.*0;
        for i = 1:2*Ns
            [~, state] = ode113(@(t, x) propagator(t, x, planetParams(1)),...
                time, [xhat_i_prev(:, i)], options);
            Xhat_i(:, i) = state(end, 1:Ns)';
        end
    
        % apriori state estimate @ time k
        Xhat_min = 1/(2*Ns).* sum(Xhat_i, 2);
    
        % apriori state covariance @ time k
        A  = zeros(Ns, Ns);
        for i = 1:2*Ns
            A = A + (Xhat_i(:, i) - Xhat_min) * (Xhat_i(:, i) - Xhat_min)';
        end
         P_min = 1/(2*Ns).* A;
    
         % create sigma points. State @ k
        [Xhat_i] = sigmaPoints_state(Ns, reshape(P_min, [36, 1]), Xhat_min);
    
        % compute predicted measurements. State @ time k
        Yhat_i = zeros(Nm, 2*Ns);
        for i = 1:2*Ns
         [Yhat_i(:, i)] = compute_measuremts(t(k), Xhat_i(:, i)', planetParams(1), zeros(6, 1));
        end
        Yhat = 1/(2*Ns).* sum(Yhat_i, 2);
    
        % compute measurement covariance
         B  = zeros(Nm, Nm);
        for i = 1:2*Ns
            B = B + (Yhat_i(:, i) - Yhat) * (Yhat_i(:, i) - Yhat)';
        end
         Py = 1/(2*Ns).* B + R0;
         
        C  = zeros(Ns, Nm);
        for i = 1:2*Ns
            C = C + (Xhat_i(:, i) - Xhat_min) * (Yhat_i(:, i) - Yhat)';
        end
         Pxy = 1/(2*Ns).* C;
    
         % kalman update
         K = Pxy/(Py);
         Xhat_plus = Xhat_min;
         P_plus = P_min - K * Py * K';
        
         % save states
         X(:, k) = Xhat_plus;
         P(k, :) = reshape(P_plus, [36, 1]);
        
         p = sqrt(diag(P_plus));
         CRB_pos(k) = max(p(1:3));
         CRB_vel(k) = max(p(4:6));
    end
end

% CR3BP propagator in the inertial frame
function [dx] = propagator(t, x, mu)
        M = t;
        r1 = [x(1)+mu*cos(M);x(2)+mu*sin(M);x(3)];              % SC-Earth
        r2 = [x(1)+(mu - 1)*cos(M); x(2)+(mu-1)*sin(M); x(3)];  % SC-Moon
        
        Cmat = zeros(7, 7); Smat = Cmat;
        Cmat(1, 1) = 1; Smat(1, 1) = 0;
        n_max = 0; normalized = 1; Re = 1;

        [~, dU1, ~] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r1, Re, 1-mu, ...
                                                    normalized);

        [~, dU2, ~] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r2, Re, mu, ...
                                                    normalized);
        dU = (dU1 + dU2);

        dx = [x(4);
          x(5);
          x(6);
          dU(1);
          dU(2);
          dU(3)];
end

% GG measurements c
function [y] = compute_measuremts(t, x, mu, noise)
        M = t;
        r1 = [x(1)+mu*cos(M);x(2)+mu*sin(M);x(3)];              % SC-Earth
        r2 = [x(1)+(mu - 1)*cos(M); x(2)+(mu-1)*sin(M); x(3)];  % SC-Moon
        
        Cmat = zeros(3, 3); Smat = Cmat;
        Cmat(1, 1) = 1; Smat(1, 1) = 0;
        n_max = 0; normalized = 1; Re = 1;

        [~, ~, ddU1] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r1, Re, 1-mu, ...
                                                    normalized);

        [~, ~, ddU2] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                    r2, Re, mu, ...
                                                    normalized);
        T = ddU1 + ddU2;
        y0 = [T(1,1); T(1, 2); T(1, 3); T(2, 2); T(2, 3); T(3, 3)];
        y = y0 + noise;
end

function [xhat_i] = sigmaPoints_state(Ns, P, Xhat)
    xhat_i = zeros(Ns, 2*Ns);
    xhat_i(:, :) = Xhat.* ones(Ns, 2*Ns);
    % compute matrix square root
    L = chol(Ns.*reshape(P, [6,6]), 'lower');
    X = L';

    % propagate first set of points
    for i = 1:Ns
        xhat_i(:, i)    =  xhat_i(:, i)    + X(i, :)';
        xhat_i(:, i+Ns) =  xhat_i(:, i+Ns) - X(i, :)';
    end
end
