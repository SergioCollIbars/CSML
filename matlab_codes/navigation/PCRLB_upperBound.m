clear;
clc;
close all;
addpath("data/")
addpath("functions/")

%%                        PCRLB UPPER BOUND
% Description: Test the sequential PCRLB in different set of trajectories
% download from JPL website and compare with the real PCRLB for different
% measurement types.
% Author: Sergio Coll Ibars
% Date: 11/13/224
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalization units
GM_E = 3.986004418E14;              % [m^3/s^2]
GM_M = 4.9048695E12;                % [m^3/s^2]
D    = 384399E3;                    % [m]
n    = sqrt((GM_E + GM_M) / D^3);   % [1/s]
mu   = GM_M / (GM_E + GM_M);        % [-]

% load initial condition and trajectory
data = load('EM_Vertical_L2_Family.mat');
X0   = data.trajFam{2,1}.iState;
P    = data.trajFam{2, 1}.period;
Ns   = 2;

% Simulation parameters
tmin = 0;           % [-]
tmax = 1*P;      % [-]
meas = "DSN";       % QGG / DSN

% integrate trajectory
options = odeset('RelTol',1e-20,'AbsTol',1e-20);
planetParams = [mu, 0, 0, 0]; poleParams = [0, 0, 0, 0];
Cmat = 0; Smat = 0; system = "CR3BP_inertial";
STM0 = reshape(eye(6,6), [36, 1]);
[t, state] = ode113(@(t, x) EOM_3BP(t, x, planetParams, ...
    poleParams, Cmat, Smat, system), [tmin, tmax], [X0; STM0], options);

% plot trajectory
figure()
plot3(state(:, 1), state(:, 2), state(:, 3), ...
    'LineWidth', 2);
axis equal;

% initial uncetainty and meas weight
sG    = 1E-12 / (n^2);                        % [-]
R_QGG = diag([sG, sG, sG, sG, sG, sG].^2);    % [-]

sR    = 1/D;                                  % [-]
sRR   = 1E-3/(D*n);                           % [-]
R_DSN = diag([sR, sRR].^2);                   % [-]

sP = 1E6/D;                                   % [-]
sV = 10/(D*n);                                % [-]
P0 = diag([sP, sP, sP, sV, sV, sV].^2);       % [-]

% compute PCRLB and upper bound to compare
P0 = P0(1:Ns, 1:Ns);
g0 = inv(P0); g = 0;

Nt  = length(t);
STM = state(:, 7:end);
pos = state(:, 1:3)';
vel = state(:, 4:6)';
PCRLB = zeros(Nt, Ns*Ns); PCRLB(1, :) = reshape(inv(P0), [Ns*Ns, 1]);
UB    = zeros(Nt, 1); CB = UB;
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
        R0 = R_QGG;
        posRel = pos(:, k) - [-mu;0;0];
        h1 = compute_grad_posPartials(1-mu, posRel(1), posRel(2), posRel(3));

        posRel = pos(:, k) - [1-mu;0;0];
        h2 = compute_grad_posPartials(mu, posRel(1), posRel(2), posRel(3));

        h = h1 + h2;
    elseif(meas == "DSN")
        R0 = R_DSN;
        posRel = pos(:, k) - [-mu;0;0];
        velRel = vel(:, k);
        h = computePartials_DSN(posRel, velRel);
        h = h(1:2, 1:Ns); 
    end

    % IF filter
    Ai_min = PHI_inv' * Aiprev_plus * PHI_inv;
    Ai_plus = Ai_min + h' * (R0 \ h);

    % get upper bopund and cov bound
    B = h' * (R0 \ h);
    g = g + det(B)^(1/Ns);
    CB(k) = (1/(det(Ai_plus)^(1/Ns)))^(Ns/2);  % [m^3]
    UB(k) = (1/(det(g0)^(1/Ns) + g)^(Ns/2));   % [m^3]
    
% %     UB(k) = ((1/min(eig(B))))^(1/2);
% %     CB(k) = max(eig(inv(Ai_plus)))^(1/2); % max direction constraint [m]

    % update information matrix
    PCRLB(k, :)=  reshape(Ai_plus, [Ns*Ns, 1]);
end

% PLOT
figure()
time = t/n/86400;                   % [days]
scale = D^3;                        % [m^3]
semilogy(time, CB.*scale , time, UB.*scale, 'LineWidth', 2)
title('Max uncertianty along orbit')
legend('truth', 'Upper bound')

% FUNCTIONS
function [h] = computePartials_DSN(posRel, velRel)
    Ns = length(posRel);
    rho_par    = [(posRel./vecnorm(posRel))', zeros(1, Ns)];
    
    a = 1/vecnorm(posRel) * velRel' * (eye(Ns,Ns) -...
        (posRel*posRel')./(vecnorm(posRel)^2));
    b = (posRel./vecnorm(posRel))';
    rhoDot_par = [a, b];
    
    h = [rho_par;rhoDot_par];
end

