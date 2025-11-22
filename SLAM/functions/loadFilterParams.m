function [R0, P0, P0c, Q0, Qb, delta_state0, Cnm, Snm, smoothing] = ...
    loadFilterParams(planetParams, instrumentParams, Cnm_t, Snm_t)
    
    % load gravity data
    loadGrav = 1;                   % options: 1 / 0

    % apply smoothing
    smoothing = 1;                  % options: 1 / 0;

    % Sampling frequency
    fs = instrumentParams(:, 5);

    % Measurement variance
    sigma_m = instrumentParams(:, 2);
    
    % Weight matrix
    R0 = diag((sigma_m.*sqrt(fs)).^2);
    
    % Process noise
    sigma_PN = [1;1;1].*1E-12;  % [m/s^2]
    Q0 = diag(sigma_PN.^2);
    
    % bias RW
    sigma_RW = 0; % [mE]
    Qb = diag((ones(6, 1).*sigma_RW).^2);

    % Filter uncertainty
    sigmaP = .1;                % [m]
    sigmaV = .01;               % [m/s]
    sigmaB = 0.5;                  % [mE]

    P0 = diag([sigmaP;sigmaP;sigmaP;sigmaV;sigmaV;sigmaV;...
        ones(6, 1).*sigmaB].^2);

    % Initial errors
    deltaP = normrnd(0, sigmaP, [3,1]); deltaV = normrnd(0, sigmaV, [3, 1]);
    delta_state0 = [deltaP; deltaV;...
        zeros(6,1)];

    % coefficient perturbation
    n_max = planetParams(3);
    [Nc, Ns, ~]  = count_num_coeff(n_max); 
    C_S_coeffs   = mat2list(Cnm_t, Snm_t, Nc, Ns);

    K = 0.025;
    n = 2:length(Cnm_t)-1;
    sigma_n = K./(n.^2);
    sigma_n(1) = sigma_n(1)*1/20;
    sigma_n(2) = sigma_n(2)*1/10;
    if(loadGrav == 0)
        [Xp, Pc]  = perturb_coeff(sigma_n, n_max, C_S_coeffs);
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, Xp);
    
        P0c  = Pc; P0c(1,1) = 0; 
    elseif(loadGrav == 1)
        Xc       = load("data/gravField_coeff.mat").Xg;
        P0c_grav = load("data/gravField_cov.mat").Pg;

        [Cnm, Snm] = list2mat(n_max, Nc, Ns, Xc);
        P0c = P0c_grav;
    end
end

