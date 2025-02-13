% % clear;
% % clc;
close all;
format short;

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

% load SPICE kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

%%                    QGG NAVIGATION FILTER
% Description: Use gradiometer measurements to position S/C in cislunar
% space. Estimation done with different filters and propagation from
% Ephemerides model
% Author: Sergio Coll Ibars
% Date: 08/20/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial configuration
system = "EPHEM";    % options: CR3BP, FCR3BP, EPHEM
solver = "EKF";      % options: CKF, EKF, batch
plotResults   = 0;   % options: 1 or 0
consider_cov  = 0;   % options: 1 or 0
process_noise = "SNC"; % options: SNC, DMC, 0
augmented_st  = 1;   % options: 1 or 0

% time parameters
tmin = 4*1.4968;
tmax = 3*1.4968 + tmin;                                 % [rad] 
frec = 1/30;                                            % measurement freq. [Hz]

% load universe
[planetParams, poleParams, Cmat_true, Smat_true, TIME, ~] = ...
    load_universe(system, [tmin, tmax], frec);
f_time = 1/30;                                         % fixed integ. frec [Hz]
n = round((TIME(end)-TIME(1))*(f_time/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);
N = f_time / frec;
DOM = ones(1, N) * NaN;
for h = 1:n/N
    val = N * (h - 1) + 1;
    DOM(h) = TIME(val);
end

% load initial conditions
<<<<<<< HEAD
% % X0_true = load_initCond(system, planetParams, TIME);
X0_true = load('initState_truth_SRP.mat').s;
=======
% % X0 = load_initCond(system, planetParams, TIME);
X0 = load('initState_truth_SRP.mat').s;
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074
STM0 = reshape(eye(6,6), [36, 1]);
if(augmented_st), X0_true = [X0_true; planetParams(10)]; STM0 = reshape(eye(7,7), [49, 1]);  end

% integrate trajectory
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[time, ~] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat_true, Smat_true, system, 0, {0,0}, augmented_st), TIME, ...
    [X0_true; STM0], options);
TIME = unique(cat(2, time', DOM));
[t, state] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
poleParams, Cmat_true, Smat_true, system, 0, {0,0}, augmented_st), TIME, ...
[X0_true; STM0], options);

% compute primaries position
[posE, posM, posS] = compute_posPrimaries(TIME, planetParams, system);

% compute gradiometer measurements
% % noiseSeed = load("noiseSeed_f10_T5.mat"). noise;
noiseSeed = [];
[T, ~, ~] = compute_measurements(TIME, state, planetParams, poleParams, ...
    Cmat_true, Smat_true, 1, 0, augmented_st, noiseSeed, DOM, posE, posM, posS, system);

% plot measurements
if(plotResults), plot_measurements(TIME, T, planetParams, augmented_st, system); end

% initialize filter. Add error to true state
[~, ~, R0, Q0, Bw, c, Pc, Pxc, Cmat_estim, Smat_estim] = ...
    initialize_filter(planetParams, Cmat_true, Smat_true, ...
    consider_cov, process_noise, augmented_st);
<<<<<<< HEAD
% % X0 = X0_true + Xnot;
=======
% % X0 = X0 + Xnot;
>>>>>>> e12bb3a6b89fc140530fa18d13ab934e4bcc0074

if(process_noise == "DMC")
    Xnot  = zeros(9, 1); X0 = [X0; zeros(3, 1)]; Ns = 9;
end
if(augmented_st == 1)
    Xnot  = zeros(7, 1); Ns = 7;
else
    Xnot  = zeros(6, 1); Ns = 6; 
end

X0 = load('initState_SRP.mat').s;
P0 = reshape(load('initCov_SRP.mat', 'p').p, [Ns,Ns]);

% compute initial error
scale = [planetParams(2)*ones(3, 1)/1E3; ...
    planetParams(2)*planetParams(3)*ones(3,1)];
err0 = (state(1, 1:6)' - X0(1:6)).*scale;   % [Km] and [m/s]

% LOOP
count = 0;                                          % iteration counter
error = 100;                                        % initial error
epsilon = 1E-15;
MaxIter = 5;
lastwarn('') % Clear last warning message 
prefIter = ones(MaxIter*6, length(t)) * NaN;
posIter  = ones(MaxIter*6, length(t)) * NaN;
corr_iter  = ones(2, MaxIter) * NaN;
error_iter = ones(2, MaxIter) * NaN;
tic
while(abs(error) > epsilon && count < MaxIter)
    % filter solver
    if(solver == "batch")
        [X, P, Xhat, XNOT, pref, posf] = batch_solver(t, X0, Xnot, P0, ...
                    R0, c, Pc, Pxc, T, planetParams, poleParams, count, ...
                    Cmat_estim, Smat_estim, system, consider_cov, augmented_st, DOM, posE, posM, posS);
        
        error_iter(1, count + 1) = vecnorm(state(1, 1:3)' - X(1:3, 1));
        error_iter(2, count + 1) = vecnorm(state(1, 4:6)' - X(4:6, 1));
        [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
                pref, posf);
    elseif(solver == "CKF")
        [X, P, Xhat, XNOT, pref, posf] = CKF_solver(t, X0, Xnot, P0, ...
                    R0, Q0, Bw, T, planetParams, poleParams, ...
                    Cmat_estim, Smat_estim, system, consider_cov, augmented_st, DOM, posE, posM, posS);

        [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
                pref, posf);
    elseif(solver == "EKF" || solver == "UKF")
        Ntmax = round(1*86400*frec);
        t_batch = t(1:Ntmax);
        while(abs(error) > epsilon && count < MaxIter) % first run batch
            [X_B, P_B, Xhat_B, XNOT, pref, posf] = CKF_solver(t_batch, X0, Xnot, P0, ...
                    R0, Q0, Bw,T, planetParams, poleParams, ...
                    Cmat_estim, Smat_estim, system, consider_cov, augmented_st, DOM, posE, posM, posS);
            
            [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
                pref, posf);
        end
        disp('Final CKF error ' + string(error))
        
        if(solver == "EKF")
            % run EKF
            t_EKF = t(1:end);
            X0 = X_B(:, 1);
            P0 = reshape(P_B(1, :), [Ns,Ns]);
            [X_E, P_E, Xhat_E, XNOT, pref, posf] = EKF_solver(t_EKF, X0, P0, ...
                        R0, Q0, Bw, T, planetParams, poleParams, ...
                        Cmat_estim, Smat_estim, system, consider_cov, augmented_st, DOM, ...
                        posE, posM, posS);
    
            X = X_E;
            P = P_E;
            Xhat = Xhat_E;
            P(1, :) = reshape(P0, [1, Ns*Ns]);

            [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
                pref, posf);
        else
            % run UKF
            t_UKF = t(1:end);
            X0 = X_B(:, 1);
            P0 = reshape(P_B(1, :), [Ns,Ns]);
            [X_U, P_U, Xhat_U, XNOT, pref, posf] = UKF_solver(t_UKF, X0, P0, ...
                        R0, Q0, Bw, T, planetParams, poleParams, ...
                        Cmat_estim, Smat_estim, system, consider_cov, augmented_st, DOM, ...
                        posE, posM, posS);
    
            X = X_U;
            P = P_U;
            Xhat = Xhat_U;
            P(1, :) = reshape(P0, [1, Ns*Ns]);

            [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
                pref, posf);
        end
    else
        Xnot = epsilon - 1;
        Xhat = ones(6, length(TIME)) * NaN;
        P = ones(length(TIME), 36) * NaN;
    end

    % catch warning if any
    [warnMsg, warnId] = lastwarn;
    if ~isempty(warnMsg)
       error = epsilon;
    end
    lastwarn('') % Clear last warning message
end
toc
disp('iterations = ' + string(count));

% clear kernels
cspice_kclear

% plot results
if(plotResults)
    plot_results(t, state, X, P, Pc, Xhat, prefIter, posIter, planetParams, ...
        count, consider_cov, augmented_st, posE, posM, system)

    % plot correction per iter
    figure()
    semilogy(1:count, corr_iter(:, 1:count), 'LineWidth', 2, 'Marker', 'sq', ...
        'LineStyle', 'none');
    title('state correction per iteration')
    xlabel('Iteration number')
    ylabel('|XNOT| [-]')
    grid on;

% %     % save data
% %     saveData(t, state, X, T);
end

function [Xnot, error, e, count, prefIter, posfIter] = ...
    check_err_save_post(Xnot, XNOT, e, count, prefIter, posfIter, pref, posf)
    % compute error
    error = vecnorm(XNOT);
    disp(error);
    e(1, count + 1) = vecnorm(XNOT(1:3));
    e(2, count + 1) = vecnorm(XNOT(4:6));
    
    if error > 1E2
        warning('Filter divergence')
    end

    % update state deviation
    Xnot = Xnot + XNOT;
    
    % Increment counter
    count = count + 1;

    % save prefit and postfit
    Nt = length(pref(1, :));
    Nm = length(pref(:, 1));
    maxInd = Nm*count;
    minInd = maxInd - (Nm-1);
    prefIter(minInd:maxInd, 1:length(pref)) = pref;
    posfIter(minInd:maxInd, 1:length(posf)) = posf;
end
