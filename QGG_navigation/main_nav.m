clear;
clc;
close all;
format short;

addpath("data/")
addpath("functions/")
addpath("simplified_functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

% load SPICE kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

%%                      GG NAVIGATION FILTER
% Description: Simplified filter code for Gravity Gradiometer (GG)
% navigation. Improves speed performance with less adaptability.
% Author: Sergio Coll Ibars
% Date: 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial configuration
system = "CR3BP";                                       % options: CR3BP, EPHEM
plotResults   = 1;                                      % options: 1 or 0
process_noise = 1;                                      % options: 1 or 0
augmented_st  = 0;                                      % options: 1 or 0

% time parameters
tmin = 0*1.4968;
tmax = 1*1.4968 + tmin;                                 % [rad] 
frec = 1/30;                                            % meas. freq. [Hz]

% load universe & initial conditions
[planetParams, poleParams, Cmat_true, Smat_true, TIME, ~] = ...
    load_universe(system, [tmin, tmax], frec);
X0_true = load_initCond(system, planetParams, TIME);

% integrate true trajectory
Ns = 14;                                                % state number:
                                                        % [r, v, eta, b]
STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
if(system == "CR3BP")
    [time, ~] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, Cmat_true, Smat_true, system, 0, {0,0}, augmented_st, 0), TIME, ...
        [X0_true; STM0], options);
elseif(system == "EPHEM")
    [time, ~] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, Cmat_true, Smat_true, system, 0, {0,0}, augmented_st, 0), TIME, ...
        [X0_true; STM0], options);
end

% compute gradiometer measurements + attitude estimates
noiseSeed = [];
[posE, posM, posS] = compute_posPrimaries(TIME, planetParams, system);
[T, ~, ~] = compute_measurements(TIME, state, planetParams, poleParams, ...
    Cmat_true, Smat_true, 1, 0, augmented_st, noiseSeed, DOM, posE, posM, posS, system);

if(plotResults), plot_measurements(TIME, T, planetParams, augmented_st, system); end

% initialize filter. Add error to a priori states
[Xnot, P0, R0, Q0, Bw, c, Pc, Pxc, Cmat_estim, Smat_estim] = ...
    initialize_filter(planetParams, Cmat_true, Smat_true, ...
    consider_cov, process_noise, augmented_st);
X0 = X0_true + Xnot;

% Estimation process
count = 0;                                          % iteration counter
error = 100;                                        % initial error
epsilon = 1E-15;
MaxIter = 5;
while(abs(error) > epsilon && count < MaxIter)
        Ntmax = round(1*86400*frec);
        t_batch = t(1:Ntmax);
        while(abs(error) > epsilon && count < MaxIter) % first run batch
            [X_B, P_B, Xhat_B, XNOT, pref, posf] = CKF_solver(t_batch, X0, Xnot, P0, ...
                    R0, Q0, Bw,T, planetParams, poleParams, Cmat_estim, Smat_estim, ...
                    system, consider_cov, augmented_st, DOM, posE, posM, posS, tau, []);
            
            [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
                pref, posf);
        end
        disp('Final CKF error ' + string(error))

        % run EKF
        t_EKF = t(1:end);
        X0 = X_B(:, 1);
        P0 = reshape(P_B(1, :), [Ns,Ns]);
        [X_E, P_E, Xhat_E, XNOT, pref, posf] = EKF_solver(t_EKF, X0, P0, ...
                    R0, Q0, Bw, T, planetParams, poleParams, ...
                    Cmat_estim, Smat_estim, system, consider_cov, augmented_st, DOM, ...
                    posE, posM, posS, bias);

        X = X_E;
        P = P_E;
        Xhat = Xhat_E;
        P(1, :) = reshape(P0, [1, Ns*Ns]);

        [Xnot, error, corr_iter, count, prefIter, posIter] = ...
            check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter, posIter, ...
            pref, posf);
       
        Xnot = epsilon - 1;
        Xhat = ones(6, length(TIME)) * NaN;
        P = ones(length(TIME), 36) * NaN;
end

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
end