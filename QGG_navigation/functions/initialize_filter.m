function [Xnot, P0, R0, Q0, Bw, c, Pc, Pxc, Cmat_estim, Smat_estim] = ...
    initialize_filter(planetParams, Cmat, Smat, consider_cov, process_noise, augmented_st, frec)
    %%                    INITIALIZE FILTER FUNCTION
    % Description: Compute inital deviation and apriori covariance for a 
    % given dynamical system.
    % Author: Sergio Coll Ibars
    % Date: 03/27/2024
    
    % apriori covariance
    % % sigmaState = [1E7, 1E7, 1E7, 10, 10, 10, 1];        % [m], [m/s] and [-]
    sigmaState = [1E3, 1E3, 1E3, 0.1, 0.1, 0.1, 1, ...
        ones(1, 6)*1E-9, ones(1, 3)*1E-3];                  % [m], [m/s],[-], [s^-2] and [m/s^2]
% %     sigmaState = [6, 6, 6, 0.001, 0.001, 0.001, 1, ...
% %         ones(1, 6)*1E-9, ones(1, 3)*1E-3];                  % [m], [m/s],[-], [s^-2] and [m/s^2]

    % apriori measurement error matrix
    sigmaMeas = [1, 1] * 1E-12 * sqrt(frec);                % [1/s^2]
    sigmaMeasAcc = ones(1, 3)  * 1E-10 * sqrt(frec);        % [m/s^2]

    % process noise. SNC
    sigmaQ_SNC  = [5E-11, 5E-11, 5E-11];                    % [m/s^2]            for PN
    sigmaQ_GG   = ones(1, 6) * 1E-16;                       % [1/s^2 * sqrt(Hz)]
    sigmaQ_acc  = ones(1, 3) * 1E-15;                       % [m/s^2 * sqrt(Hz)] 

    sigma_Bw    = 30 / (6*86400);                           % STD BW noise [1/s] for Jupiter
    sigmaQ_DMC  = [1E-5, 1E-5, 1E-5];                       % [m/s^2]            for Jupiter

    % inital deviation
% %     XNOT = [1E6; 1E6; 1E6; 5; 5; 5; 0.5];               % [m], [m/s] and [-]
    XNOT = [1E1; 1E1; 1E1; .05; .05; .05; 0.5;...
        ones(6, 1)*1E-10; ones(3, 1)*1E-4];                 % [m], [m/s],[-],[s^-2] and [m/s^2]
% %     XNOT = [1; 1; 1; 1E-4; 1E-4; 1E-4; 0.5;...
% %         ones(6, 1)*1E-12; ones(3, 1)*1E-6];                 % [m], [m/s],[-],[s^-2] and [m/s^2]


    % consider parameters & consider uncertainty
    path_sc1 = "SIGMACOEFS_EARTH_1.txt";
    path_sc2 = "SIGMACOEFS_MOON_1.txt";
    sc1 = readmatrix(path_sc1);
    sc2 = readmatrix(path_sc2);

    n_max = planetParams(6);
    n_data = n_max;
    [Nc, Ns] = countCoeff(n_max);
    [Nc_data, Ns_data] = countCoeff(n_data);

    sc1 = sc1(4:end).*1E0;
    sc2 = sc2(4:end).*1E0;
    [X1] = mat2list(Cmat{1}, Smat{1}, Nc_data, Ns_data);
    [X2] = mat2list(Cmat{2}, Smat{2}, Nc_data, Ns_data);
    X1  = [X1(1:Nc); X1(Nc_data+1:Ns+Nc_data)];
    X2  = [X2(1:Nc); X2(Nc_data+1:Ns+Nc_data)];

% %     % generate estimation coeffiecients
% %     sigma1  = [sc1(1:Nc); sc1(Nc_data+1:Ns+Nc_data)];
% %     sigma2  = [sc2(1:Nc); sc2(Nc_data+1:Ns+Nc_data)];
% %     s1_rand = normrnd(0, sigma1(2:end));
% %     s2_rand = normrnd(0, sigma2(2:end));

    % load consider parameters random error seed
    s1_rand = load("considerSeed_1_n8.mat").s1_rand;
    s2_rand = load("considerSeed_2_n8.mat").s2_rand;
    
    if(consider_cov)
        X1c = X1;
        X2c = X2;
        X1c(2:end) = X1(2:end) + s1_rand;
        X2c(2:end) = X2(2:end) + s2_rand;

        % coefficient errors
        c = [X1c(2:end)-X1(2:end);X2c(2:end)-X2(2:end)].*0;

        Pc = diag([sc1(2:Nc);sc1(Nc_data+1:Ns+Nc_data);...
          sc2(2:Nc);sc2(Nc_data+1:Ns+Nc_data)].^2);
    else
        X1c = X1;
        X2c = X2;
        X1c(2:end) = X1(2:end) + s1_rand;
        X2c(2:end) = X2(2:end) + s2_rand;

        % coefficient errors
        c  = 0;
        Pc = 0;
    end
    [Cmat1, Smat1] = list2mat(X1c, n_max);
    [Cmat2, Smat2] = list2mat(X2c, n_max);
    Cmat_estim = {Cmat1, Cmat2};
    Smat_estim = {Smat1, Smat2};
   
    Pxc = zeros(6, length(c));

    % transform to non-dimensional units
    distDim = planetParams(2);
    velDim  = planetParams(2) * planetParams(3);
    measDim = planetParams(3)^2;
    measDimAcc = planetParams(2) * planetParams(3)^2;
    Hz      = planetParams(3);

    sigmaState(1:3)  = sigmaState(1:3)./distDim;     % [-]
    sigmaState(4:6)  = sigmaState(4:6)./velDim;      % [-]
    sigmaState(8:13) = sigmaState(8:13)./measDim;    % [-]
    sigmaState(14)   = sigmaState(14)./measDimAcc;   % [-]
    if(augmented_st)
        P0 = diag(sigmaState.^2);                    % [-]
    else
        P0 = diag(sigmaState(1:6).^2);               % [-]
    end
    
    if(process_noise == "SNC") 
        sigmaQ_SNC = sigmaQ_SNC./(measDim*distDim);                   % [-]
        Q0 = diag([sigmaQ_SNC(1), sigmaQ_SNC(2), sigmaQ_SNC(3)].^2);  % [-]
        % % Bw = zeros(3,3);

        sigmaQ_GG  = sigmaQ_GG./(planetParams(3)^2 * sqrt(Hz));
        sigmaQ_acc = sigmaQ_acc./(measDimAcc * sqrt(Hz));
        Q_GG       = diag(sigmaQ_GG.^2);
        Q_acc      = diag(sigmaQ_acc.^2);
        Bw         = [Q_GG, zeros(6,3);zeros(3,6),Q_acc];
    elseif(process_noise == "DMC")
        sigma_Bw = sigma_Bw/planetParams(3);                        % [-]
        sigmaQ_DMC = sigmaQ_DMC/(measDim*distDim);                  % [-]

        QTilde0 = diag([sigmaQ_DMC(1), sigmaQ_DMC(2), sigmaQ_DMC(3)].^2);   % [-]
        Bw = eye(3, 3) * sigma_Bw;                                          % [-]
        P0 = [P0, zeros(6, 3);zeros(3, 6), QTilde0];                        % [-]
        Q0 = QTilde0;
    else
        Q0 = zeros(3,3);
        Bw = zeros(3,3);
    end

    if(augmented_st)
        sigmaMeas = sigmaMeas./(measDim);                % [-]
        sigmaMeasAcc = sigmaMeasAcc./(measDimAcc);       % [-]
        R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
        sigmaMeas(2), sigmaMeas(1), sigmaMeasAcc(1), sigmaMeasAcc(2), ...
        sigmaMeasAcc(3)].^2);                            % [-]

        Xnot(1:3, 1) = XNOT(1:3)./distDim;               % [-]
        Xnot(4:6, 1) = XNOT(4:6)./velDim;                % [-]
        Xnot(7, 1)   = XNOT(7);                          % [-]
        Xnot(8:13, 1)= XNOT(8:13)./measDim;              % [-]
        Xnot(14:16, 1)  = XNOT(14:16)./measDimAcc;       % [-]
    else
        sigmaMeas = sigmaMeas./(measDim);               % [-]
        R0 = diag([sigmaMeas(1), sigmaMeas(2), sigmaMeas(2), sigmaMeas(1), ...
        sigmaMeas(2), sigmaMeas(1)].^2);                % [-]

        Xnot(1:3, 1) = XNOT(1:3)./distDim;              % [-]
        Xnot(4:6, 1) = XNOT(4:6)./velDim;               % [-]
    end
end