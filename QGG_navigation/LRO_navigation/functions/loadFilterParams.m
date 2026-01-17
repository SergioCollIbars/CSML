function [R0, P0, Q0, Qb, delta_state0, Cnm, Snm, instrument_alig, sigma_att] = ...
    loadFilterParams(folder_Name, planetParams, instrumentParams, ...
    Cnm_list, Snm_list)
    p = readParams("data/"+folder_Name+"/Filter.txt");
    
    % Instrument aligment
    instrument_alig = p.orientation;

    % Sampling frequency
    fs = instrumentParams(:, 5);

    % Measurement std
    sigma_m = [p.sigma_xx;p.sigma_xy;p.sigma_xz;...
             p.sigma_yy;p.sigma_yz;p.sigma_zz];
    
    % Weight matrix
    R0 = diag((sigma_m.*sqrt(fs)).^2);
    
    % Process noise
    sigma_PN = [p.sigma_PN_x;p.sigma_PN_y;p.sigma_PN_z];  % [m/s^2]
    Q0 = diag(sigma_PN.^2);
    
    % bias RW
    sigma_RW = p.sigma_RW;               % [mE]
    Qb = diag((ones(6, 1).*sigma_RW).^2);

    % attitude uncertainty
    sigma_att = p.sigma_att * pi / (180 * 3600);           % [radians]

    % Filter uncertainty
    sigmaP = p.sigmaP;                   % [m]
    sigmaV = p.sigmaV;                   % [m/s]
    sigmaB = p.sigmaB;                   % [mE]

    P0 = diag([sigmaP;sigmaP;sigmaP;sigmaV;sigmaV;sigmaV;...
        ones(6, 1).*sigmaB].^2);

    % Initial errors
    deltaP = normrnd(0, sigmaP, [3,1]); 
    deltaV = normrnd(0, sigmaV, [3, 1]);
    deltaB = normrnd(0, sigmaB, [6, 1]);
    delta_state0 = [deltaP; deltaV;...
        deltaB];

    % coefficient perturbation (Earth & Moon)
    Cnm_E  = Cnm_list{1}; Snm_E = Snm_list{1};
    Cnm_M  = Cnm_list{2}; Snm_M = Snm_list{2};
    n_max        = planetParams(5);
    [Nc, Ns, ~]  = count_num_coeff(n_max); 
 
    C_S_coeffs_E   = mat2list(Cnm_E, Snm_E, Nc, Ns);
    C_S_coeffs_M   = mat2list(Cnm_M, Snm_M, Nc, Ns);

    K = 0;
    n_E = 2:length(Cnm_E)-1; n_M = 2:length(Cnm_M)-1;
    sigma_n_E = K./(n_E.^2); sigma_n_M = K./(n_M.^2);

    [Xp, ~]  = perturb_coeff(sigma_n_E, n_max, C_S_coeffs_E);
    [Cnm_E, Snm_E] = list2mat(n_max, Nc, Ns, Xp);

    [Xp, ~]  = perturb_coeff(sigma_n_M, n_max, C_S_coeffs_M);
    [Cnm_M, Snm_M] = list2mat(n_max, Nc, Ns, Xp);
    
    % save coefficients as a list 
    Cnm = {Cnm_E, Cnm_M}; Snm = {Snm_E, Snm_M};
end

