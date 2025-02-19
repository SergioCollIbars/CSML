clear;
clc;
close all;
addpath("data/")
addpath("functions/")
set(0,'defaultAxesFontSize',16);

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

% load initial condition and trajectory
data = load('EM_Lyap_L2_Family.mat');
% % data = load('EM_DRO_Family.mat');
% % data = load('EM_LoPO_Family.mat');
% % data = load('EM_31_Res_Family.mat');
Nd =length(data.trajFam);
index = 1:40:Nd;

% trajectories struct
dataOrbit = struct('traj', nan(3, 1E4), 'time', nan(1, 1E4), 'JC', []);
dataCRB_pos    = struct('value', nan(1, 1E4));
for i = 1:length(index)
    % Filling with NaN data 
    dataOrbit(i).traj = ones(3, 1E4)*NaN;
    dataOrbit(i).time = ones(1, 1E4)*NaN;
    dataOrbit(i).JC   = data.trajFam{i, 1}.CJ; 
    dataCRB_pos(i).value   = ones(1, 1E4) * NaN;
end
dataCRB_vel = dataCRB_pos;

for j = 1:length(index)
    x0   = data.trajFam{index(j),1}.iState;
    mu = data.trajFam{index(j),1}.mu;
    P    = data.trajFam{index(j), 1}.period;
    Ns   = 3;
    
    [r0, v0] = rotate2inertial(x0(1:3), x0(4:6), 0, 1);
    X0 = [r0;v0];
    
    % Simulation parameters
    tmin = 0;           % [-]
    tmax = 1*P;         % [-]
    frec = 1/30;        % [Hz]
    TIME = linspace(tmin, tmax, round(tmax*frec/n));
    meas = "QGG";       % QGG / DSN
    
    % integrate trajectory
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);
    planetParams = [mu, 0, 0, 0]; poleParams = [0, 0, 0, 0];
    Cmat = 0; Smat = 0; system = "CR3BP_inertial";
    STM0 = reshape(eye(6,6), [36, 1]);
    [t, state] = ode113(@(t, x) EOM_3BP(t, x, planetParams, ...
        poleParams, Cmat, Smat, system), TIME, [X0; STM0], options);
    dataOrbit(j).traj = state(:, 1:3)';
    dataOrbit(j).time = t';
    
    % compute Earth and Moon positions in inertial frame
    [posE, posM] = compute_EM_position(t, mu);
    
    % initial uncetainty and meas weight
    sG    = 1E-12 / (n^2);                         % [-]
    R_QGG = diag([sG, sG, sG, sG, sG, sG].^2);     % [-]
    
    sP = 1E20/D;                                   % [-]
    sV = 100/(D*n);                                % [-]
    P0 = diag([sP, sP, sP, sV, sV, sV].^2);        % [-]
    
    % compute Upper bounds
    scaleP = D;     % [m]
    scaleV = D*n;   % [m/s]
    [CRB_pos, CRB_vel] = compute_CRB(state, posM, posE, t, P0, R_QGG, Ns,...
        mu, meas);
    dataCRB_pos(j).value = CRB_pos.*scaleP;    % [m]
    dataCRB_vel(j).value = CRB_vel.*scaleV;    % [m/s]
end

% plot trajectory
figure()
colormap("jet");
maxVal = ones(1, length(index)) * NaN;
minVal = maxVal;
for j = 1:length(index)
[rB] = rotate2BodyFrame(dataOrbit(j).time,  dataOrbit(j).traj');
scalarValues = dataCRB_pos(j).value;
maxVal(j) = max(scalarValues);
minVal(j) = min(scalarValues);
scatter3(rB(1, :), rB(2, :), rB(3, :), 20, log10(scalarValues), 'filled');
axis equal;
hold on;
end
plot(-mu, 0, "o",'MarkerFaceColor',"#7E2F8E", 'MarkerEdgeColor', "#7E2F8E")
plot((1-mu),0, "o",'MarkerFaceColor',"#77AC30", 'MarkerEdgeColor', "#77AC30")
colorbar; % Show color scale
% % caxis([log10(min(minVal)) log10(max(maxVal))]); % Adjust color range to scalar values
% % caxis([-3, 6])
title('Maximum position value from covariance [m]')

% plot trajectory
figure()
colormap("jet");
maxVal = ones(1, length(index)) * NaN;
minVal = maxVal;
for j = 1:length(index)
[rB] = rotate2BodyFrame(dataOrbit(j).time,  dataOrbit(j).traj');
scalarValues = dataCRB_vel(j).value;
maxVal(j) = max(scalarValues);
minVal(j) = min(scalarValues);
scatter3(rB(1, :), rB(2, :), rB(3, :), 20, log10(scalarValues), 'filled');
axis equal;
hold on;
end
plot(-mu, 0, "o",'MarkerFaceColor',"#7E2F8E", 'MarkerEdgeColor', "#7E2F8E")
plot((1-mu),0, "o",'MarkerFaceColor',"#77AC30", 'MarkerEdgeColor', "#77AC30")
colorbar; % Show color scale
% % caxis([log10(min(minVal)) log10(max(maxVal))]); % Adjust color range to scalar values
% % caxis([-3, 6])
title('Maximum velocity value from covariance [m/s]')

figure()
for j =1:length(index)
    semilogy(dataOrbit(j).time, dataCRB_pos(j).value, 'LineWidth', 2)
    ylabel('[m]')
    hold on;
end

% FUNCTIONS
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

function [CRB_pos, CRB_vel] = compute_CRB(state, posM, posE, t, P0, R0, Ns, mu, meas)
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
        PHI = PHI_2 * inv(PHI_1);               % from k-1 to k
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
             CRB_vel(k)  = 1;                  % max direction constraint   [m/s]
        end
    
        % update information matrix
        PCRLB(k, :)=  reshape(Ai_plus, [Ns*Ns, 1]);
    end
end

