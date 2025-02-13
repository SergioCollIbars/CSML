clear;
clc;
close all;
format long g;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
%%                    QGG NAVIGATION OBSERVABILITY
% Description: compute the information matrix and upper bpund at different
% locations in space.
% Author: Sergio Coll Ibars
% Date: 11/12/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% BOUNDARY TEST
clear;
clc;
close all;
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

% specify simulation time, frequency and noise
tmin = 0;                                          % Phi = 0  [-]
% % tmin = 0.75449775462963;                           % Phi = 180[-]
tmax = 3.4968 + tmin;                              % [-]
frec = 1/60;

% define system
system = "CR3BP"; % options: 2BP, CR3BP, F2BP
[planetParams, ~, ~, ~, ~, ~] = ...
    load_universe(system, [tmin, tmax], frec);

% normalization values
measDim = planetParams(3)^2;
timeDim = planetParams(3);
posDim  = planetParams(2);                            % [m]
velDim  = planetParams(3) * planetParams(2);          % [m/s]

% Measurement weights
sigmaMeas = [1, 1/2] * 1E-12;                        % [1/s^2]
sigmaMeas = sigmaMeas./measDim;                      % [-]
sa = 1E-9 / (planetParams(2) * planetParams(3)^2);   % [-]
R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
sigmaMeas(2), sigmaMeas(1)].^2);                     % [-]

% load parameters & initial conditions
[planetParams, poleParams, C_mat, S_mat, TIME, ~] = ...
    load_universe(system, [tmin, tmax], frec);
X0 = load_initCond(system, planetParams);

[~, P0, ~, ~, ~, ~, ~, ~, ~, ~] = ...
    initialize_filter(planetParams, C_mat, S_mat, ...
    0, "0", 0);

% compute measurement time
f_time = frec;
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);
N = f_time / frec;
DOM = ones(1, N) * NaN;

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6,6), [36, 1]);
[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, C_mat, S_mat, system, 0, {0,0}, 0), TIME, [X0; STM0], options);

STM  = state(:, 7:42);
Nt  = length(TIME);

% primaries motion
[posE, posM] = compute_posPrimaries(TIME, planetParams, system);

% compute PCRB
A = ones(1, Nt)*NaN; B = A;
Ns = 3; n = Ns;
PCRB = ones(Nt, Ns*Ns) * NaN; P = PCRB;
PCRB(1, :) = reshape(inv(P0(1:3, 1:3)), [1, Ns*Ns]);
g0 = inv(P0(1:3, 1:3));
g = 0;
for k = 1:Nt
    if(k == 1)
        PHI_1 = eye(6,6);
        Aiprev_plus = reshape(PCRB(k, :), [Ns,Ns]);
    else
        PHI_1 = reshape(STM(k-1, :), [6,6]);            % from 0 to k-1
        Aiprev_plus = reshape(PCRB(k-1, :), [Ns,Ns]);
    end
    PHI_2 = reshape(STM(k, :), [6,6]);                  % from 0 to k
    PHI = PHI_2* inv(PHI_1);       % from k-1 to k;
    PHI = PHI(1:3, 1:3);

    % compute measurements and visibility matrix for gradiometer
    [~, Hx, ~] = compute_measurements(TIME(k), state(k, 1:6), planetParams, ...
         poleParams, C_mat, S_mat, 0, 0, 0, [], DOM, posE(:, k), posM(:, k), [], system);

    h = Hx(:, 1:3);

    % gamma function
    PHI_inv = inv(PHI);
    
    M = PHI_inv' * Aiprev_plus * PHI_inv;
    Ai_min = M;
    
    Ai_plus = Ai_min + h' * (R0 \ h);
    f = h' * inv(R0) * h;

    alpha = 1;
    l = min(eig(f));

    % max direction constraint
    B(k) = ((1/alpha) * (1/l))^(1/2);
    A(k) = max(eig(inv(Ai_plus)))^(1/2); % max direction constraint

% %     % volume constraint
% %     B(k) = 1/alpha^(3/2) * (1/det(f))^(1/2);
% %     B(k) = (1/(det(g)^(1/3) + det(f)^(1/3)))^(1/2);
% %     A(k) = (1/det(Ai_plus))^(1/6);           
    
    g = g + det(f)^(1/n);
    A(k) = (1/(det(Ai_plus)^(1/n)))^(n/2);  % [m^3]
    B(k) = (1/(det(g0)^(1/n) + g)^(n/2));   % [m^3]
    
    % store new value
    PCRB(k, :) = reshape(Ai_plus, [1, Ns*Ns]);
end

figure()
scale = planetParams(2);    % [m]
% % scale = 1;
semilogy(t./timeDim/86400, A.*scale, ...
    t./timeDim/86400, B.*scale, 'LineWidth', 2)
legend('truth cov.', 'upper bounding')
xlabel('Time [days]')
ylabel('[m]')

% clear kernels
cspice_kclear
%%              COMPUTE INFO. AT DIFFERENT POINTS

% Domain definiton. Non-dimensional units 
Nx = 200;
Ny = 200;
IF_QGG   = ones(Nx, Ny) * NaN;
BD_QGG   = ones(Nx, Ny) * NaN;

% % minD = -0.05;
% % maxD = 0.05;
% % [X, Y] = meshgrid(linspace(minD, maxD, Nx), linspace(-0.2, maxD, Ny));

% % minD = -0.01;
% % maxD = 0.01;
% % [X, Y] = meshgrid(linspace(minD, maxD, Nx), linspace(minD, maxD, Ny));

minD = -1.5;
maxD = 1.5;
% % [X, Y] = meshgrid(linspace(minD, maxD, Nx), linspace(minD, maxD, Ny));
[X, Y] = meshgrid(linspace(-1.1, 1.3, Nx), linspace(-1, 1, Ny));

% measurements matrices
sigmaMeas = [1, 1/2] * 1E-12;                           % [1/s^2]
sigmaMeas = sigmaMeas./(planetParams(3)^2);             % [-]
R_QGG = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
    sigmaMeas(2), sigmaMeas(1)].^2);                    % [-]
s = 1E-9 /(planetParams(2) * planetParams(3)^2);        % [-]
R_acc = diag([s, s, s].^2);                             % [-]
mu = planetParams(1);

for j = 1:Nx    % x domain
    for i =1:Ny % y domain
        time = 0;   % [-]
        state = [X(j, i), Y(j, i), 0];
% %         state = [1-mu, X(j, i), Y(j, i)];

        % gradiometer sensitivity
        [posE, posM] = compute_posPrimaries(time, planetParams, system);

        [~, Hc, ~] = compute_measurements(time, state, planetParams, poleParams, ...
            C_mat, S_mat, 0, 0, 0, [], DOM, posE, posM, [], system);

        % visibility matrix
        h = Hc(:, 1:3);
        R = R_QGG;

        % rescale measurements
        scale  = planetParams(2);   % [m]

        f = h' * inv(R) * h;
        l = min(eig(f));
        BD = ((1/l))^(1/2);                           % bounding for largest uncertanty direction
        IF  = det(f);
        BD_QGG(j, i)   = BD * scale;                  % [m]
        IF_QGG(j, i)   = IF * 1/scale^3;              % [1/m^3]
    end
end

scale  = planetParams(2)/1E3;   % [km]
lg = ["", "",  "Earth", "Moon", "L_1", "L_2", "L_3", "L_4", "L_5"];
ttitle = "Information gain value (\delta |\Delta_x|) in log10 scale";
plot_sensitivityContour(X, Y, log10(IF_QGG), ttitle, ...
    planetParams(1), scale, lg);

ttitle = "PCRLB upper bound around Moon. Gradiometer measurement";
plot_sensitivityContour(X, Y, log10(BD_QGG), ttitle, ...
    planetParams(1), scale, lg);


% clear kernels
cspice_kclear

%% FUNCTIONS
function [] = plot_sensitivityContour(X, Y, Z, ttitle, mu, scale, lg)
    % Lagrange points
    L_1 = .836915 * scale;
    L_2 = 1.15568 * scale;
    L_3 = -1.00506 * scale;
    
    % trajectory
    figure()
    contourf(X.*scale, Y.*scale, Z, 'EdgeColor', 'none');
    view(0, 90);
    colorbar()
    xlabel('X [Km]')
    ylabel('Y [Km]')
    hold all;
    axis equal;
    plot(-mu*scale,0, "o",'MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
    plot((1-mu)*scale,0,"pentagram",'MarkerFaceColor','r', 'MarkerEdgeColor', 'r')
    plot(L_1,0,'rv','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot(L_2,0,'r^','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot(L_3,0,'rp','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot((0.5-mu)*scale,(sqrt(3)/2)*scale,'rX','MarkerFaceColor','k', 'MarkerEdgeColor', 'k')
    plot((0.5-mu)*scale,(-sqrt(3)/2)*scale,'rs','MarkerFaceColor','k', 'MarkerEdgeColor','k')
    legend('Earth', 'Moon', 'L1', 'L2', 'L3', 'L4', 'L5')
    legend(lg)
    title(ttitle)
end