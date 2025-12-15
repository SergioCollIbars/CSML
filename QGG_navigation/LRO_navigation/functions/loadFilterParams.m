function [R0, P0, P0c, Q0, Qb, delta_state0, Cnm, Snm, Nbatch] = ...
    loadFilterParams(metaData_file, planetParams, instrumentParams, ...
    Cnm_t, Snm_t)
    mtd = readParams("data/"+metaData_file);
    p = readParams("data/"+mtd.folder+"/Filter.txt");
    
    % load gravity data
    loadGrav = p.loadGrav;

    % Batch arc for gravit field estimation [sec]
    t_batch = p.t_batch * 3600;           % [sec]

    % Sampling frequency
    fs = instrumentParams(:, 5);

    % Batch arc smaples number
    Nbatch = round(t_batch * fs(1));

    % Measurement variance
    sigma_m = instrumentParams(:, 2);
    
    % Weight matrix
    R0 = diag((sigma_m.*sqrt(fs)).^2);
    
    % Process noise
    sigma_PN = [p.sigma_PN_x;p.sigma_PN_y;p.sigma_PN_z];  % [m/s^2]
    Q0 = diag(sigma_PN.^2);
    
    % bias RW
    sigma_RW = p.sigma_RW;               % [mE]
    Qb = diag((ones(6, 1).*sigma_RW).^2);

    % Filter uncertainty
    sigmaP = p.sigmaP/3;                   % [m]
    sigmaV = p.sigmaV/3;                   % [m/s]
    sigmaB = p.sigmaB;                     % [mE]

    P0 = diag([sigmaP;sigmaP;sigmaP;sigmaV;sigmaV;sigmaV;...
        ones(6, 1).*sigmaB].^2);

    % Initial errors
    deltaP = normrnd(0, sigmaP, [3,1]); 
    deltaV = normrnd(0, sigmaV, [3, 1]);
    deltaB = normrnd(0, sigmaB, [6, 1]);
    delta_state0 = [deltaP; deltaV;...
        deltaB];

    % coefficient perturbation
    n_max        = planetParams(3);
    [Nc, Ns, ~]  = count_num_coeff(n_max); 
    C_S_coeffs   = mat2list(Cnm_t, Snm_t, Nc, Ns);

    K = 0.043;
% %     K = 1E-20;
    n = 2:length(Cnm_t)-1;
    sigma_n = K./(n.^2);
    if(loadGrav == 0)
        [Xp, Pc]  = perturb_coeff(sigma_n, n_max, C_S_coeffs);
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, Xp);
    
        P0c  = Pc; P0c(1,1) = 0; 
    elseif(loadGrav == 1)
        Xc       = load("data/"+mtd.gravity).Xg;
        P0c_grav = load("data/"+mtd.cov).Pg;

        [Cnm, Snm] = list2mat(n_max, Nc, Ns, Xc);
        P0c = (1.3^2).*P0c_grav;
    end
end

