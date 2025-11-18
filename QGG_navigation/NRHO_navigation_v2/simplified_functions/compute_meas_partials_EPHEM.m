function [Hmeas] = compute_meas_partials_EPHEM(t, state, planetParams,...
    BN, C_mat, S_mat, posE, posM, posS)
    % Compute measurements partials using the EPHEMERIDES model. 
    % The measurements partials computed are: gradiometer partials
    % and acceloremeter partials.
    % Date: 09/24/2025
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % S/C mass parameters
    mass = planetParams(11);    % [Kg]
    A    = planetParams(12);    % [-]
    eta  = 1.3;
    
    % pointer to current time index
    j = 1;

    % gradiometer partials
    [H_GG_pos] = compute_posPartials(planetParams, BN, C_mat, S_mat, t(j),...
                state(j, 1:3)', posE(:, j), posM(:, j), posS(:, j));
    H_GG_pos   = H_GG_pos * (planetParams(3)^2);      % [s^-2 / [-]]

    H_GG_vel   = zeros(6, 3);
    H_GG_eta   = zeros(6, 1);
    H_GG_bias  = [eye(6, 6)];                         % [E]
    
    r3         = state(j, 1:3)' - posS(:, j);
    
    [~, daSRP_dr, daSRP_deta] = SRP(r3, eta, mass, A);
    daSRP_dr = BN * daSRP_dr * BN';
    daSRP_deta = BN * daSRP_deta;
            
    % accelerometer partials
    H_acc_pos  = daSRP_dr * 1E3;
    H_acc_vel  = zeros(3, 3);
    H_acc_eta  = daSRP_deta * 1E3;
    H_acc_bias = [zeros(3, 6)]; 
    
    % ensamble measurement matrix
    measDim = 1E-9;
    Hmeas = [H_GG_pos/measDim*1E-3, H_GG_vel, H_GG_bias];
end

