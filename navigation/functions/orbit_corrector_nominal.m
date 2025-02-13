function [state_true, state_nom, state_recons, state_error_true_nom, state_error_true_recons] =...
    orbit_corrector_nominal...
    (X0, delta, planetParams, poleParams ,Cmat, Smat, t, n_max)
    
    % define integration options
    options = odeset('RelTol',1e-13,'AbsTol',1e-13);

    GM = planetParams(1);
    Re = planetParams(2);
    normalized = planetParams(4);
    
    % ODE 113
    [~, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
        poleParams, Cmat, Smat), t, X0, options);
    [~, state_nom] = ode113(@(t, x) EOM_navigation(t, x, [GM, Re, n_max, normalized], ...
        poleParams, Cmat, Smat), t, X0+delta, options);
    
    % generate gradiometer measurements
    [ddU_ACI] = gradiometer_meas(t, state_true(:, 1:3), Cmat, Smat, ...
        poleParams, planetParams);
    
    % reconstruct trajectory
    delta = (state_true(2, :) - state_true(1, :))';
    [delta_new] = reconstruct_traj(t, ddU_ACI, delta);
    
    % state reconstruction from nominal
    state_recons = zeros(6, length(t));
    state_recons(:, 1) = X0;
    for j = 2:length(t)
        state_recons(:, j) = state_recons(:, j-1) + delta_new(:, j);
    end
    
    % error states
    state_error_true_nom = state_true - state_nom;
    state_error_true_recons = state_true - state_recons';
end

