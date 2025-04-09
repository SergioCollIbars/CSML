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
%%                    QGG NAVIGATION OBSERVABILITY
% Description: compute the observability of the S/C state in the desired
% system using gradiometer measurements. Observability is function of the
% meas frequency and noise.
% Author: Sergio Coll Ibars
% Date: 04/24/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
clear;
close all;
clc;

% specify simulation time, frequency and noise
tmin = 0;                                            % Phi = 0  [-]
% % tmin = 0.75449775462963;                         % Phi = 180[-]
% % tmin =  0.615368240740741;                       % Phi = 30 [-]
T = 1.4968;                                          % [-]
tmax = 1*T + tmin;                                   % [-]
frec = 1/30;
augmented_st = 0;
bias  = 1;

% define system
system = "CR3BP"; % options: 2BP, CR3BP, F2BP
[planetParams, ~, ~, ~, ~, ~] = ...
    load_universe(system, [tmin, tmax], frec);

% normalization values
measDim = planetParams(3)^2;
timeDim = planetParams(3);                            % [1/s]
posDim  = planetParams(2);                            % [m]
velDim  = planetParams(3) * planetParams(2);          % [m/s]

% load parameters & initial conditions
[planetParams, poleParams, C_mat, S_mat, TIME, ~] = ...
    load_universe(system, [tmin, tmax], frec);
X0 = load_initCond(system, planetParams);
% % X0 = [ 0.720230409024078;0.675424403969901;0.00837326048601556;...
% %     -1.81989246557446;1.92454671114969;0.119861341700249];  % state @ periapsis
% % X0 = [0.838678033240026; 0.542367087602379 ; -0.0740369975323576;...
% %     -0.624211548421406;0.826244659528416;0.42190501223657];     % state @ phi = 30

% compute measurement time
f_time = 1/30;
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);
N = f_time / frec;
DOM = ones(1, N) * NaN;
for h = 1:n/N
    val = N * (h - 1) + 1;
    DOM(h) = TIME(val);
end
dt = TIME(2) - TIME(1);                 % [-]

% Measurement weights
sigmaMeas = [1, 1/sqrt(2)] * 1E-12;                  % [1/s^2]
sigmaMeas = sigmaMeas./measDim;                      % [-]
sAcc      = 1E-10./(planetParams(2)*planetParams(3)^2); % [-] 
R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
sigmaMeas(2), sigmaMeas(1), sAcc, sAcc, sAcc].^2);    % [-]

%  state process noise
sigmaQ_s = 5E-8/ (planetParams(2)*planetParams(3)^2) ;% [-]
qs = diag([sigmaQ_s, sigmaQ_s, sigmaQ_s].^2).*1;
I = eye(3, 3);
Nq = 6;

% random-walk bias process noise (FOGM)
sigmaQ = 1E-5;                                       % [E / sec]
varQ   = sigmaQ^2;                                   % [E^2 / sec^2]
qE     = varQ * dt / timeDim;                        % [E^2 / sec] assuming it is multiplied by dt
q      = qE * 1E-18;                                 % [1 / sec^5]
q      = q  / (timeDim^5);                           % [-]
tau    = 5 * timeDim;                                % [-]

if(augmented_st), Nx = 7; Ns = 7; Nm = 6; 
elseif(bias), Nx = 18; Ns = 6; Nm = 6;
else, Nx = 6;  Ns = 6; Nm = 6; end

if(tau == 0)
    Nx  = 12; % constant bias
    q = 0;
elseif(tau/timeDim >= 1E20 )
    Nx  = 18; % random-walk bias
    PHIbd = [1, dt;0, 1];
    Sbd   = zeros(2,2);
else
    PHIbd = compute_STM_FOGM(dt, tau);         
    Sbd   = compute_PN_FOGM(dt, tau, q);   
end

if(augmented_st), X0 = [X0; planetParams(10)]; end
P0 = diag(repelem(1E3, Nx));

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, C_mat, S_mat, system, 0, {0,0}, augmented_st), TIME, ...
    [X0; STM0], options);

STM  = state(:, Ns+1:Ns+Ns*Ns);
Nt  = length(TIME);

% primaries motion
[posE, posM, posS] = compute_posPrimaries(TIME, planetParams, system);

% compute PCRB
PCRB = ones(Nt, Nx*Nx) * NaN; P = PCRB;
sigmaP = ones(Nx, Nt) * NaN;
PCRB(1, :) = reshape(inv(P0), [1, Nx*Nx]);
obs = ones(1, Nt) * NaN;
for k = 2:Nt
    PHI_1 = reshape(STM(k-1, :), [Ns, Ns]); % from 0 to k-1
    PHI_2 = reshape(STM(k, :), [Ns, Ns]);   % from 0 to k
    PHI = PHI_2 * inv(PHI_1);               % from k-1 to k;

    % compute measurements and visibility matrix
    [Y, Ht, ~] = compute_measurements(TIME(k), state(k, 1:Ns), planetParams, ...
         poleParams, C_mat, S_mat, 0, 0, augmented_st, [], DOM, posE(:, k), ...
         posM(:, k), posS(:, k), system);
    
    % gamma function
     At = t(k) - t(k-1);
     Gamma = At * [At/2*I;I];
     PHI_inv = inv(PHI);
     Aiprev_plus = reshape(PCRB(k-1, :), [Nx,Nx]);
     if(Nx == 18) % include bias noise     
         F = diag_stack(PHI, PHIbd, PHIbd, PHIbd, PHIbd, PHIbd, PHIbd);
         Qy = diag_stack(Gamma * qs * Gamma', Sbd, Sbd, Sbd, Sbd, Sbd, Sbd);
         M = inv(F)' * Aiprev_plus * inv(F);
         Ai_min = M - M*Qy*inv(eye(Nx, Nx) + M * Qy)*M;
         Hb = compute_biasDrift_measPartials();
         h = [Ht, Hb];
     else % constant bias or random-walk       
         F = [PHI, zeros(6,6);zeros(6,6), eye(6,6)];
         Qq = q * dt * eye(6,6);
         Qy = diag_stack(Gamma * qs * Gamma', Qq);
         M = inv(F)' * Aiprev_plus * inv(F);
         Ai_min = M - M*Qy*inv(eye(Nx, Nx) + M * Qy)*M;
         
         h = [Ht, eye(6, 6)];
     end
    if(isnan(vecnorm(Y)))
        Ai_plus = Ai_min;
    else
        Ai_plus = Ai_min + h' * (R0(1:Nm, 1:Nm) \ h);
    end
    obs(k) = rank(Ai_plus);

    % store new value
    PCRB(k, :) = reshape(Ai_plus, [1, Nx*Nx]);
end

% compute & propagate state uncertanty
for j = 1:Nt
    p = inv(reshape(PCRB(j, :), [Nx,Nx])); 
    sigmaP(:, j)   = sqrt(diag(p));
    sigmaP(1:3, j) = sigmaP(1:3, j).*(posDim/1E3);
    sigmaP(4:6, j) =  sigmaP(4:6 ,j).*(velDim);
end

% plot observability
figure()
plot(TIME/timeDim/86400, obs, 'LineWidth', 2)
xlabel('TIME [days]')
ylabel('[-]');
title('system observability')

% plot uncertainty
figure()
sfig = [1, 3, 5, 2, 4, 6];
tt = ["X position", "Y position", "Z position", "X velocity", "Y velocity", ...
    "Z velocity"];
tty = ["[Km]", "[Km]", "[Km]", "[m/s]", "[m/s]", "[m/s]"];
for j = 1:6
    subplot(3, 2, sfig(j))
    semilogy(TIME/timeDim/86400, 3*sigmaP(j, :), 'LineWidth', 2)
    title(tt(j));
    xlabel('TIME [days]')
    ylabel(tty(j))
end
tt = ["b_{xx}", "b_{xy}", "b_{xz}", "b_{yy}", "b_{yz}", ...
    "b_{zz}"];
tty = ["[E]", "[E]", "[E]", "[E]", "[E]", "[E]"];
if(Nx == 18)
    b = zeros(6, length(TIME)); d = b;
    for j = 1:Nt
        b(:, j) = [sigmaP(7, j); sigmaP(9, j); sigmaP(11, j); sigmaP(13, j);...
            sigmaP(15, j); sigmaP(17, j)];
        d(:, j) = [sigmaP(8, j); sigmaP(10, j); sigmaP(12, j); sigmaP(14, j);...
            sigmaP(16, j); sigmaP(18, j)];
    end
    figure()
    subplot(1, 2, 1)
    semilogy(TIME/timeDim/86400, b.*measDim/1E-9, 'LineWidth', 2)
    title('Bias uncertianty [Eotvos]')
    subplot(1, 2, 2)
    semilogy(TIME/timeDim/86400, d.*timeDim^3/1E-9, 'LineWidth', 2)
    title('Drift uncertianty [Eotvos / sec]')
   
else
     b = zeros(6, length(TIME)); d = b;
    for j = 1:Nt
        b(:, j) = [sigmaP(7, j); sigmaP(8, j); sigmaP(9, j); sigmaP(10, j);...
            sigmaP(11, j); sigmaP(12, j)];
    end
    figure()
    semilogy(TIME/timeDim/86400, b.*measDim/1E-9, 'LineWidth', 2)
    title('Bias uncertianty [Eotvos]')
end
% clear kernels
cspice_kclear

%%                  FUNCTIONS

function [HOmega_pos, HOmega_vel] = compute_angVel_partials(r, v)    
    h = cross(r, v);

    r_skew = compute_skewMatrix(r);
    v_skew = compute_skewMatrix(v);

    HOmega_vel = +1/(vecnorm(r)^2) * r_skew;
    HOmega_pos = -1/(vecnorm(r)^2) * v_skew - 2/(vecnorm(r)^4) * ...
        [r(1)*h(1), r(2)*h(1), r(3)*h(1);...
         r(1)*h(2), r(2)*h(2), r(3)*h(2);...
         r(1)*h(3), r(2)*h(3), r(3)*h(3)];
end

function [x_skew] = compute_skewMatrix(x)
    x_skew  = [0, -x(3), x(2);...
               x(3), 0, -x(1);...
               -x(2), x(1), 0];
end

function [HOmegaDot_pos, HOmegaDot_vel] = compute_angAcc_partials(r, v, dU, ddU)
    h = cross(r, v);
    hdot = cross(r, dU);
    rn = vecnorm(r);
    rh = r./rn;
    rn2 = rn^2;
    rn3 = rn^3;
    rn4 = rn^4;
    rn5 = rn^5;

    r_skew = compute_skewMatrix(r);
    v_skew = compute_skewMatrix(v);
    a_skew = compute_skewMatrix(dU);

    HOmegaDot_vel = -2/rn3 * (v'*rh) * r_skew - 2/rn4 * ...
        [r(1)*h(1), r(2)*h(1), r(3)*h(1);...
         r(1)*h(2), r(2)*h(2), r(3)*h(2);...
         r(1)*h(3), r(2)*h(3), r(3)*h(3)];

    HOmegaDot_pos = -a_skew/rn2 - 2/rn4 * ...
        [r(1)*hdot(1), r(2)*hdot(1), r(3)*hdot(1);...
         r(1)*hdot(2), r(2)*hdot(2), r(3)*hdot(2);...
         r(1)*hdot(3), r(2)*hdot(3), r(3)*hdot(3)] -2/rn4 * ...
        [v(1)*h(1), v(2)*h(1), v(3)*h(1);...
         v(1)*h(2), v(2)*h(2), v(3)*h(2);...
         v(1)*h(3), v(2)*h(3), v(3)*h(3)] + 2/rn3 * (v' * rh) * v_skew + ...
         6/rn5 * (v' * rh) * [r(1)*h(1), r(2)*h(1), r(3)*h(1);...
         r(1)*h(2), r(2)*h(2), r(3)*h(2);...
         r(1)*h(3), r(2)*h(3), r(3)*h(3)] + 1/rn2 * r_skew * ddU;
end

function [Htot] = arrangePartials(Hgrad, Homega, HomegaDot, w)
    Htot = ones(6, 3);

    % xx direction partials
    Htot(1, :) = Hgrad(1, :) -2*w(2)*Homega(2, :) - 2*w(3)*Homega(3, :);
    % xy direction partials
    Htot(2, :) = Hgrad(2, :) + w(2)*Homega(1, :) + w(1)*Homega(2, :) - HomegaDot(3, :);
    % xz direction partials
    Htot(3, :) = Hgrad(3, :) + w(3)*Homega(1, :) + w(1)*Homega(3, :) + HomegaDot(2, :);
    % yy direction partials
    Htot(4, :) = Hgrad(4, :) -2*w(3)*Homega(3, :) - 2*w(1)*Homega(1, :);
    % yz direction partials
    Htot(5, :) = Hgrad(5, :) + w(3)*Homega(2, :) + w(2)*Homega(3, :) - HomegaDot(1, :);
    % zz direction partials
    Htot(6, :) = Hgrad(6, :) -2*w(1)*Homega(1, :) - 2*w(2)*Homega(2, :);
end


function BigMatrix = diag_stack(varargin)
    % Determine total size
    total_size = sum(cellfun(@(M) size(M, 1), varargin));
    
    % Initialize the big matrix
    BigMatrix = zeros(total_size);
    
    % Fill the diagonal blocks
    idx = 1;
    for i = 1:length(varargin)
        M = varargin{i};
        s = size(M, 1);
        BigMatrix(idx:idx+s-1, idx:idx+s-1) = M;
        idx = idx + s;
    end
end

function [Hb] = compute_biasDrift_measPartials()
    Hb = zeros(6, 12);
    for j = 1:6
        maxInd = j * 2;
        minInd = maxInd - 1;
        Hb(j, minInd:maxInd) = [1, 0];
    end
end

function [PHIbd] = compute_STM_FOGM(t, tau)
    PHIbd = [1, tau*(1 - exp(-t/ tau));...
        0, exp(-t/tau)];
end

function [Sbd] = compute_PN_FOGM(t, tau, sigmaQ_b)
    S11 = tau^2 * ((1-exp(-2*t/tau)) - 4*(1-exp(-t/tau)) + 2*t/tau);
    S12 = tau * (1 - exp(-t/tau))^2;
    S21 = tau * (1 - exp(-t/tau))^2;
    S22 = (1 - exp(-2*t/tau));

    Sbd = sigmaQ_b * tau / 2 * [S11,S12;S21,S22];
end