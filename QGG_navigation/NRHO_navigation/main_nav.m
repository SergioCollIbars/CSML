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
% navigation. Improves speed performance with less adaptability. Works with
% Epehemrides model only.
% Author: Sergio Coll Ibars
% Date: 09/24/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initial configuration
plotResults   = 1;                                      % options: 1 or 0
saveData      = 0;                                      % options: 1 or 0
loadData      = 0;                                      % options: 1 or 0
attitude      = "inertial";                             % options: inertial

% time parameters
tmin = 0;
tmax = 2*1.4968 + tmin;                                 % [rad] 
frec = 1/30;                                            % meas. freq. [Hz]

if(loadData), [tmin,~,~,~] = loadData_files(); tmax = tmax + tmin; end

% load universe & initial conditions
[planetParams, poleParams, Cmat_true, Smat_true, TIME, ~] = ...
    load_universe("EPHEM", [tmin, tmax], frec);

n = round((TIME(end)-TIME(1))*(frec/planetParams(3)) + 1);
TIME = linspace(TIME(1), TIME(end), n);

% % TIME = TIME /planetParams(3);
% % planetParams(2) = 1; planetParams(3) =1;

traj0 = load_initCond("EPHEM", planetParams, TIME);
[b0_GG, b0_acc] = load_bias(planetParams);

if(loadData == 1)
    [~,X0_true,~,~] = loadData_files(); 
elseif(loadData == 0)
    X0_true = [traj0; planetParams(10);b0_GG;b0_acc]; 
end

% integrate true trajectory
disp('Intregrating true trajectory ...')
Ns = 16;                                                % state number:
                                                        % [r, v, eta, b]
STM0 = reshape(eye(Ns,Ns), [Ns*Ns, 1]);
options = odeset('RelTol',1e-12,'AbsTol',1e-12);
[time, state] = ode113(@(t, x) EOM_navigation_EPHEM(t, x, planetParams, ...
    Cmat_true, Smat_true), TIME, [X0_true; STM0], options);
disp('  DONE!')

% compute gradiometer measurements + attitude estimates
disp('Computing measurements ...')
[posE, posM, posS] = compute_posPrimaries(TIME, planetParams, "EPHEM");
[BN_matrix_true] = compute_attitude_EPHEM(TIME, frec, 1, "inertial");
[BN_matrix_meas] = compute_attitude_EPHEM(TIME, frec, 0, "inertial");
[T, acc] = compute_measurement_EPHEM(TIME, state, planetParams, ...
    BN_matrix_true, Cmat_true, Smat_true, 1, posE, posM, posS);

if(plotResults), plot_measurements(TIME, [T;acc], ...
        planetParams, 1, "EPHEM"); end
disp('  DONE!')

% initialize filter. Add error to a priori states
[Xnot, P0, R0, Q0, Qb, ~, ~, ~, Cmat_estim, Smat_estim] = ...
    initialize_filter(planetParams, Cmat_true, Smat_true, ...
    0, "SNC", 1, frec);
if(loadData == 1)
    [~,~,X0, P0] = loadData_files(); 
elseif(loadData == 0)
    X0 = X0_true + Xnot;
end
Xnot = Xnot.*0;

% Estimation process
count = 0;                                          % iteration counter
error = 100;                                        % initial error
epsilon = 1E-15;
MaxIter = 5;

prefIter = ones(MaxIter*6, length(TIME)) * NaN;
posIter  = ones(MaxIter*6, length(TIME)) * NaN;
corr_iter  = ones(2, MaxIter) * NaN;
error_iter = ones(2, MaxIter) * NaN;
while(abs(error) > epsilon && count < MaxIter)
        if(loadData == 0)
            Ntmax = round(1*86400*frec);
% %             Ntmax = round(length(TIME) * 0.1);
            t_batch = TIME(1:Ntmax);
            while(abs(error) > epsilon && count < MaxIter) % first run CKF to initialize
                [X_B, P_B, Xhat_B, XNOT, pref, posf] = ...
                    CKF_solver_EPHEM(t_batch, X0, Xnot, P0, R0, Q0, Qb,...
                    [T;acc], planetParams, BN_matrix_meas, ...
                    Cmat_estim, Smat_estim, posE, posM, posS);
                
                [Xnot, error, corr_iter, count, prefIter, posIter] = ...
                    check_err_save_post(Xnot, XNOT, corr_iter, count,...
                    prefIter, posIter, pref, posf);
            end
            disp('Final CKF error ' + string(error))
        elseif(loadData == 1)
            X_B = X0;
            P_B = reshape(P0, [1, Ns*Ns]);
        end

        % run EKF
        t_EKF = TIME(1:end);
        X0 = X_B(:, 1);
        P0 = reshape(P_B(1, :), [Ns,Ns]);
        [X_E, P_E, Xhat_E, XNOT, pref, posf] = EKF_solver_EPHEM(t_EKF, X0, P0, ...
                    R0, Q0, Qb, [T;acc], planetParams, BN_matrix_meas, ...
                    Cmat_estim, Smat_estim, posE, posM, posS);

        X = X_E;
        P = P_E;
        Xhat = Xhat_E;
        P(1, :) = reshape(P0, [1, Ns*Ns]);

        [Xnot, error, corr_iter, count, prefIter, posIter] = ...
            check_err_save_post(Xnot, XNOT, corr_iter, count, prefIter,...
            posIter, pref, posf);
end

% clear kernels
cspice_kclear

% plot results
if(plotResults)
    plot_results_EPHEM(TIME, state, X, P, prefIter, posIter,...
        planetParams, count);

    % plot correction per iter
    figure()
    semilogy(1:count, corr_iter(:, 1:count), 'LineWidth', 2, 'Marker', 'sq', ...
        'LineStyle', 'none');
    title('state correction per iteration')
    xlabel('Iteration number')
    ylabel('|XNOT| [-]')
    grid on;
end

% save results
if(saveData == 1)
    X_true_final  = state(end, 1:16);
    X_estim_final = X(:, end);
    t_final       = tmax;
    P_final       = P(end, :);
    
    save('GG_navigation_final_params.mat',...
        't_final','X_true_final','X_estim_final','P_final');
end


%% FUNCTIONS
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
    Nm = length(pref(:, 1));
    maxInd = Nm*count;
    minInd = maxInd - (Nm-1);
    prefIter(minInd:maxInd, 1:length(pref)) = pref;
    posfIter(minInd:maxInd, 1:length(posf)) = posf;
end