clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Compare NSM vs LS.
% Author: Sergio Coll
% Date: 10/28/25

target_body = "Bennu";  % options: Bennu or Eros
mc_num      = 1;
Solver      = "LS_NSM"; % options: LS / NSM / LS_NSM / LS_BNSM_ANSM
nuisance    = "pos";    % options: pos / att / both

if(target_body == "Bennu")
    % Asteroid parameters.
    path = "HARMCOEFS_BENNU_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 5.2;
    n_max  = 6;
    normalized = 1;
    W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
    W0 = 0;                   % Initial asteroid longitude
    RA = deg2rad(86.6388);    % Right Ascension     [rad]
    DEC = deg2rad(-65.1086);  % Declination         [rad]
    r = 0.36E3;               % orbit radius        [m]
    sigmaP = 0.1;             % orbit error         [m]
    sigmaAtt = 10 * pi / (180*3600); % Euler angle error    [rad]
    sigmaW = 1E-3;                    % attitude errors     [deg/sqrt(hr)]
elseif(target_body == "Eros")
    path = "HARMCOEFS_EROS_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    n_max  = 10;
    normalized = 1;
    GM =  459604.431484721;          % Point mass value    [m^3/s^2]
    W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
    W0 = 0;                          % Initial asteroid longitude
    RA = deg2rad(11.363);            % Right Ascension     [rad]
    DEC = deg2rad(17.232);           % Declination         [rad]
    r = 24E3;                        % orbit radius        [m]
    sigmaP = 3;                      % orbit error         [m]
    sigmaAtt = 10 * pi / (180*3600); % Euler angle error [rad]
    sigmaW = 1E-3;                   % attitude errors     [deg/sqrt(hr)]
end
poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

Cnm = Cnm(1:n_max+1, 1:n_max+1);
Snm = Snm(1:n_max+1, 1:n_max+1);
% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 3;
f   = 1/10;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);
At = t(2) - t(1);

% Integrate trajectory (true trajectory)
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;STM0], options);

% perturb nominal coefficient
[X]   = mat2list(Cnm, Snm, Nc, Ns);
sigmaCoeff = 1;

% Gradiometer measurement noise
sigma  = 1E-11 * sqrt(f);
R  =  diag([1,1/2,1/2,1,1/2,1])*sigma^2;
RR =  diag([1,1/2,1/2,1/2,1,...
    1/2,1/2,1/2,1])*sigma^2;
STD = sqrt(diag(RR)');
mu = zeros(length(STD), 1);
Nm = length(RR);

noise0 = zeros(9, Nt);

dX1 = ones(Ncs, mc_num) * NaN; dX3 = ones(Ncs, mc_num) * NaN;
dX2 = ones(Ncs, mc_num) * NaN;
for Mc = 1:mc_num
    disp('Mc run = ' + string(Mc))
    
    % attitude errors
    deltaAtt = normrnd(0, sigmaAtt, [3, Nt]) + [1;6;10].*(pi/(180*3600));
    sigmaW   = sigmaW / sqrt(At) * pi / (180*3600);   % [rad/s]
    deltaW   = normrnd(0, sigmaW, [3, Nt]);

    % position errors
    deltaP  = normrnd(0, sigmaP, [3, Nt])   + [.01;.01;.01];

    % % deltaP   = [.1;.1;.1].*ones(3, Nt);
    % % deltaAtt = [1;6;10].*(pi/(180*3600)).*ones(3, Nt);

    % create gradiometer measurements
    Y_true = ones(9, length(t));
    noise = normrnd(repmat(mu, 1, Nt), repmat(STD(:), 1, Nt));
    for k = 1:length(t)
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        r = state_t(k, 1:3)'; v = state_t(k, 4:6)';
        [Y_ACI, ~] = gradiometer_meas(t(k) ,[GM, Re, n_max, normalized],...
            ACAF_ACI, [r', v'], zeros(9, 1), Cnm, Snm);
        
        % attitude state & angular motion error
        if(nuisance == "att" || nuisance == "both")
            BN = rotationMatrix(deltaAtt(1, k), deltaAtt(2, k),...
                deltaAtt(3, k), [3, 2, 1]);
            omega    = zeros(3, 1) + deltaW(:, k);
        else
            BN = eye(3);
            omega    = zeros(3, 1);
        end
         % compute angular velocity dyad
        [~, omega2_dyad]   = skewMatrix(omega);

        T_ACI = reshape(Y_ACI, [3,3]);
        T_B   = BN * T_ACI * BN' + omega2_dyad;
        Y_true(:, k) = reshape(T_B, [9, 1]);
    end
    % include noise
    Y_true = Y_true + noise;
    
    % compute nominal gravity field
    deltaX = normrnd(0, sigmaCoeff, [Ncs, 1]);
    deltaX(1) = 0;
    [Cp, Sp] = list2mat(n_max, Nc, Ns, X+deltaX);

    % generate initial nominal orbit
    rn = state_t(:, 1:3)';
    vn = state_t(:, 4:6)';
    if(nuisance == "pos" || nuisance == "both")
        rn = rn + deltaP;
    end
    
    % select consider parameters
    if(nuisance == "pos")
        Pc = eye(3).*(sigmaP)^2;
        Pxc = zeros(Ncs - 1, 3);
        c   = deltaP(:, 1);
    elseif(nuisance == "att")
        Pc = eye(3).*(sigmaAtt)^2;
        Pxc = zeros(Ncs - 1, 3);
        c   = deltaAtt(:, 1);
    elseif(nuisance == "both")
        Pc = diag([ones(3, 1).*sigmaP;ones(3, 1).*sigmaAtt].^2);
        Pxc = zeros(Ncs - 1, 6);
        c   = [deltaP(:, 1);deltaAtt(:, 1)];
    end
    
    if(Solver == "LS")
        P0_LS = eye(Ncs-1).*(100)^2;
        
        [X_LS, P_LS] = solver_LS(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_LS, Pc, Pxc, c, Y_true, R, nuisance);
            
        P3  = P_LS(1:Ncs-1,1:Ncs-1);
        dX3(:, Mc) = X_LS;
    elseif(Solver == "NSM")
        P0_NSM = eye(Ncs-1).*(100)^2;
        [X_NSM, P_NSM] = solver_NSM(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_NSM, Y_true, R, nuisance);
        
        P1  = P_NSM; 
        dX1(:, Mc) = X_NSM; 
    elseif(Solver == "LS_NSM")
        % solve LS
        disp('  Solve LS')
        P0_LS = eye(Ncs-1).*(100)^2;
        
        [X_LS, P_LS] = solver_LS(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_LS, Pc, Pxc, c, Y_true, R, nuisance);
            
        P3  = P_LS(1:Ncs-1,1:Ncs-1);
        dX3(:, Mc) = X_LS;

        % Solve NSM
        disp('  Solve NSM')
        P0_NSM = eye(Ncs-1).*(100)^2;
        [X_NSM, P_NSM] = solver_NSM(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_NSM, (0.3).*Pc, Pxc, Y_true, R, nuisance);
        
        P1  = P_NSM; 
        dX1(:, Mc) = X_NSM; 
    elseif(Solver == "LS_BNSM_ANSM")
        % Solve BNSM
        disp('  Solve BNSM')
        P0_NSM = eye(Ncs-1).*(100)^2;
        [X_NSM, P_NSM] = solver_BNSM(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_NSM, Y_true, R, nuisance);
        
        P2  = P_NSM; 
        dX2(:, Mc) = X_NSM; 

        % solve LS
        disp('  Solve LS')
        P0_LS = eye(Ncs-1).*(100)^2;
        
        [X_LS, P_LS] = solver_LS(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_LS, Pc, Pxc, c, Y_true, R, nuisance);
            
        P3  = P_LS(1:Ncs-1,1:Ncs-1);
        dX3(:, Mc) = X_LS;

        % Solve ANSM
        disp('  Solve ANSM')
        P0_NSM = eye(Ncs-1).*(100)^2;
        [X_NSM, P_NSM] = solver_ANSM(t, rn, vn, Cp, Sp, asterParams, ...
            poleParams, P0_NSM, Y_true, R, nuisance);
        
        P1  = P_NSM; 
        dX1(:, Mc) = X_NSM; 
    end
end

% Compute actual errors
Nx  = length(deltaX);

% Compute correlation matrix
if(Solver == "NSM")
    std_devs = sqrt(diag(P1));
    R1 = P1./ (std_devs * std_devs');
    
    % Plot correlations
    figure;
    imagesc(abs(R1));
    colormap jet;
    colorbar;
    axis equal tight;
    xticks(1:size(R1,1));
    yticks(1:size(R1,1));
    xlabel('Variable Index');
    ylabel('Variable Index');
    title('Correlation Matrix. NSM');

    Xerr1 = abs(dX1 - X);
    sigma1 = sqrt(diag(P1));

    figure()
    semilogy(1:Nx-1, 3.*sigma1, 'Color', 'b', 'LineWidth', 2);
    hold all;
    for Mc = 1:mc_num
        semilogy(1:Nx-1, Xerr1(2:Nx, Mc), 'Color', 'b', 'Marker','*', ...
            'LineStyle','none')
        hold on;
    end
    legend('\sigma NSM', 'err NSM');
    grid on;
elseif(Solver == "LS")
    std_devs = sqrt(diag(P3));
    R3 = P3./ (std_devs * std_devs');

    % Plot correlations
    figure;
    imagesc(abs(R3));
    colormap jet;
    colorbar;
    axis equal tight;
    xticks(1:size(R3,1));
    yticks(1:size(R3,1));
    xlabel('Variable Index');
    ylabel('Variable Index');
    title('Correlation Matrix. Co-Estimation method');

    Xerr3 = abs(dX3 - X);
    sigma3 = sqrt(diag(P3));

    figure()
    semilogy(1:Nx-1, 3.*sigma3, 'Color', 'r', 'LineWidth', 2);
    hold all;
    for Mc = 1:mc_num
        semilogy(1:Nx-1, Xerr3(2:Nx, Mc), 'Color', 'r', 'Marker','*', ...
            'LineStyle','none');
        hold on;
    end
    legend('\sigma LS', 'err LS');
    grid on;
elseif(Solver == "LS_NSM" || Solver == "LS_BNSM_ANSM" )
    Xerr1 = abs(dX1 - X);
    sigma1 = sqrt(diag(P1));

    if(Solver == "LS_BNSM_ANSM")
        Xerr2 = abs(dX2 - X);
        sigma2 = sqrt(diag(P2));
    end

    Xerr3 = abs(dX3 - X);
    sigma3 = sqrt(diag(P3));

    % order true values
    [X_order, ~]     = orderValues(abs(X(2:end)), n_max);
    
    % order sigma values
    [sigma1_ord, xvals] = orderValues(sigma1, n_max);
    [sigma3_ord, ~]     = orderValues(sigma3, n_max);

    figure()
    semilogy(xvals, X_order, 'Color', 'k','LineWidth', 2);
    hold all;
    semilogy(xvals, 3.*sigma1_ord, 'Color', 'b', 'LineWidth', 2);
    semilogy(xvals, 3.*sigma3_ord, 'Color', 'g', 'LineWidth', 2);
    if(Solver == "LS_BNSM_ANSM")
        [sigma2_ord, ~]     = orderValues(sigma2, n_max);
        semilogy(xvals, 3.*sigma2_ord, 'Color', 'r', 'LineWidth', 2);
    end
    
    alpha = 0.3;    % transparency
    for Mc = 1:mc_num
        % order error values
        [err1, ~] = orderValues(Xerr1(2:end, Mc), n_max);
        [err3, ~] = orderValues(Xerr3(2:end, Mc), n_max);

        % plot
        scatter(xvals, err1, 40, 'b', '*', ...
            'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha)
        hold all;
        scatter(xvals, err3, 40, 'g', '*', ...
            'MarkerFaceAlpha', alpha, 'MarkerEdgeAlpha', alpha);
        if(Solver == "LS_BNSM_ANSM")
            [err2, ~] = orderValues(Xerr2(2:end, Mc), n_max);
            semilogy(xvals, err2, 'Color', 'r', 'Marker','*', ...
            'LineStyle','none');
        end
    end
    semilogy(xvals, 3.*sigma1_ord, 'Color', 'b', 'LineWidth', 2);
    hold all;
    semilogy(xvals, 3.*sigma3_ord, 'Color', 'g', 'LineWidth', 2);
    if(Solver == "LS_BNSM_ANSM")
        legend('truth', '3 \sigma NSM', '3 \sigma LS', '3 \sigma BNSM');
    else
        legend('truth', '3 \sigma NSM', '3 \sigma LS');
    end
    grid on;
    xlabel('Degree')
    title('Gravity coefficients estimation for ' + string(target_body));
    xticks(2:n_max);
    xlim([2, n_max+1]);
    if(target_body == "Bennu")
        ylim([1E-10 1E-2]);
    else
        ylim([1E-8 1E-1]);
    end
end

%% FUNCTIONS
function [X0_final, P0_final] = solver_LS(t, rn, vn, Cnm, Snm, asterParams, ...
    poleParams, P0, Pc, Pxc, c, Y_true, R, nuisance)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 

    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns)];

    count = 0;
    iterMax = 5;
    err = 0;
    xnot_L = zeros(Ncs-1, 1);
    while (count < iterMax) && (err < 1E3)
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;
        [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
        for j = 1:Nt        
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt = W0 + W * t(j);
            ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
                [rn(:, j)', vn(:, j)'], ...
                    zeros(9,1), Cnm, Snm);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
                Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
        
             % compute position and attitude partials
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

            % compute Euler angles partials
            [Hatt] = compute_rotPartials_analy(Y_ACI, B_ACI);
            Hatt = [Hatt(1:3, :); Hatt(5:6, :);Hatt(9, :)];

            % compute consider parameters
            if(nuisance == "pos")
                [mxc, mcc] = compute_considerParams(R, Hc, Hpos);
            elseif(nuisance == "att")
                [mxc, mcc] = compute_considerParams(R, Hc, Hatt);
            elseif(nuisance == "both")
                [mxc, mcc] = compute_considerParams(R, Hc, [Hpos,Hatt]);
            end

            % prefit residuals
            dY = Y_true(:, j) - Y_ACI;
            dy = [dY(1:3);dY(5:6);dY(9)];

            % solve normal equations
            ax = Hc' * (R\Hc);
            nx = Hc' * (R\dy);
            
            Ax_L  = Ax_L + ax;
            Nx_L  = Nx_L + nx;
            Mxc = Mxc + mxc;
            Mcc = Mcc + mcc;
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L - Ax_L\(Mxc * c);
    
        X_new(2:end) = X_new(2:end) + XNOT_L;
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  LS update = '    + string(err));
    
        % update counter
        count = count + 1;
    end
    
    Px  = inv(Ax_L);
    Sxc = -Px * Mxc;
    Pxx = Px + Sxc*Pc*Sxc';
    Pxc = Sxc * Pc;

    X0_final = X_new;
    P0_final = [Pxx, Pxc;Pxc', Pc];
end

function [X_final, P_final] = solver_NSM(t, rn, vn, Cnm, Snm, asterParams, ...
    poleParams, P0, Pc, Pxc, Y_true, R, nuisance)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    
    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns)];

    count = 0;
    iterMax = 5;
    err = 0;
    xnot_L = zeros(Ncs-1, 1);
    while (count < iterMax) && (err < 1E3)
        % integrate state
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;
        [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
        for j = 1:Nt        
            % current position
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt = W0 + W * t(j);
            ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
                [rn(:, j)', vn(:, j)'], ...
                    zeros(9,1), Cnm, Snm);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
                Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
        
             % compute position and attitude partials
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

             % compute Euler angles partials
            [Hatt] = compute_rotPartials_analy(Y_ACI, B_ACI);
            Hatt = [Hatt(1:3, :); Hatt(5:6, :);Hatt(9, :)];

            if(nuisance == "pos")
                Hn = Hpos;
            elseif(nuisance == "att")
                Hn = Hatt;
            elseif(nuisance == "both")
                Hn = [Hpos, Hatt];
            end

            % Compute NS
            c = null(Hn');
            if(isempty(c))
                [~,v,D] = svd(Hn');
                nv = length(diag(v));
                C = D(:, nv-1);
                res = C' * Hn;
            else
                C = c;
                res = zeros(length(c(1, :)), length(c(1, :)));
            end
            
            % visibility matrix
            H = C' * Hc;
            r = C' * R * C;

            % prefit residuals
            dY = Y_true(:, j) - Y_ACI;
            dy = C' * [dY(1:3);dY(5:6);dY(9)];
            
            % compute residuals
            mxc = (H'   * (r\res));
            mcc = (res' * (r\res)); 

            % solve normal equations
            ax = H' * (r\H);
            nx = H' * (r\dy);
            
            Ax_L  = Ax_L + ax;
            Nx_L  = Nx_L + nx;

            Mxc = Mxc + mxc;
            Mcc = Mcc + mcc;
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X_new(2:end) = X_new(2:end) + XNOT_L;
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  NSM update = '    + string(err));
    
        % update counter
        count = count + 1;
    end

    X_final = X_new;
    Px = inv(Ax_L);
    Sxc = -Px * Mxc;
    Pxx = Px + Sxc*Pc*Sxc';
    Pxc = Sxc * Pc;
    P_N =  [Pxx, Pxc;Pxc', Pc];
    Ncs = length(Pxx);
    P_final = P_N(1:Ncs,1:Ncs);
end

function [X_final, P_final] = solver_BNSM(t, rn, vn, Cnm, Snm, asterParams, ...
    poleParams, P0, Y_true, R, nuisance)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    
    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns)];

    count = 0;
    iterMax = 8;
    err = 0;
    xnot_L = zeros(Ncs-1, 1);
    
    index  = 0;
    Hc = [];
    Hn = [];
    Dy = [];
    while (count < iterMax) && (err < 1E3)
        % integrate state
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;
        for j = 1:Nt        
            % current position
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt = W0 + W * t(j);
            ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
                [rn(:, j)', vn(:, j)'], ...
                    zeros(9,1), Cnm, Snm);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
                Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
        
             % compute position and attitude partials
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

             % compute Euler angles partials
            [Hatt] = compute_rotPartials_analy(Y_ACI, B_ACI);
            Hatt = [Hatt(1:3, :); Hatt(5:6, :);Hatt(9, :)];

            if(nuisance == "pos")
                hn = Hpos;
            elseif(nuisance == "att")
                hn = Hatt;
            elseif(nuisance == "both")
                hn = [Hpos, Hatt];
            end
            
            % prefit residuals
            res = Y_true(:, j) - Y_ACI;
            dy  = [res(1:3);res(5:6);res(9)];

            index = index + 1;
            if(index == 2)
                Hn = [Hn;hn];
                Dy = [Dy;dy];
                Hc = [Hc;hc];
                
                var = diag(R);
                R_comb = diag([var; var]);

                % Compute NS
                C = null(Hn');
    
                % visibility matrix
                H  = C' * Hc;
                r  = C' * R_comb * C;
                DY = C' * Dy;
    
                % solve normal equations
                ax = H' * (r\H);
                nx = H' * (r\DY);
                
                Ax_L  = Ax_L + ax;
                Nx_L  = Nx_L + nx;

                % reset matrices
                Hn = [];
                Dy = [];
                Hc  = [];
                index = 0;
            else
                Hn = [Hn;hn];
                Dy = [Dy;dy];
                Hc = [Hc;hc];
            end
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X_new(2:end) = X_new(2:end) + XNOT_L;
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  BNSM update = '    + string(err));
    
        % update counter
        count = count + 1;
    end

    X_final = X_new;
    P_final = inv(Ax_L);
end

function [X_final, P_final] = solver_ANSM(t, rn, vn, Cnm, Snm, asterParams, ...
    poleParams, P0, Y_true, R, nuisance)
    % extract asteroid parameteres
    GM = asterParams(1); Re = asterParams(2); n_max = asterParams(3);
    normalized = asterParams(4); Nt = length(t);

    % extract asteroid pole parameters
    W = poleParams(1); W0 = poleParams(2); RA = poleParams(3);
    DEC = poleParams(4);

    % coefficient numbers
    [Nc, Ns, Ncs] = count_num_coeff(n_max); 
    
    % Initial state
    X_new  = [mat2list(Cnm, Snm, Nc, Ns)];

    count = 0;
    iterMax = 5;
    err = 0;
    xnot_L = zeros(Ncs-1, 1);
    while (count < iterMax) && (err < 1E3)
        % integrate state
        Ax_L = inv(P0);
        Nx_L = -inv(P0) * xnot_L;
        for j = 1:Nt        
            % current position
            rn_ACI = rn(:, j);
    
            % ACAF to ACI rotation matrix
            Wt = W0 + W * t(j);
            ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
            B_ACI     = eye(3,3);
            ACAF_BODY = ACAF_ACI * B_ACI';
    
            % computed meas & select measurements
            [Y_ACI, Hc_ACI] = gradiometer_meas(t(j) ,asterParams, ACAF_ACI,...
                [rn(:, j)', vn(:, j)'], ...
                    zeros(9,1), Cnm, Snm);
            Hc_BODY = rotate_coeffPartials(Hc_ACI, B_ACI);
            Hc = [Hc_BODY(1, 2:end); Hc_BODY(4, 2:end); Hc_BODY(7, 2:end); ...
                Hc_BODY(5, 2:end);...
                Hc_BODY(8, 2:end); Hc_BODY(9, 2:end)];
        
             % compute position and attitude partials
            [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);
            Hpos = [Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)];

             % compute Euler angles partials
            [Hatt] = compute_rotPartials_analy(Y_ACI, B_ACI);
            Hatt = [Hatt(1:3, :); Hatt(5:6, :);Hatt(9, :)];

            if(nuisance == "pos")
                Hn = Hpos;
            elseif(nuisance == "att")
                Hn = Hatt;
            elseif(nuisance == "both")
                Hn = [Hpos, Hatt];
            end

            % Compute NS
            [~,v,D] = svd(Hn');
            nv = length(diag(v));
            C = D(:, nv-1);

            % visibility matrix
            H = C' * Hc;
            r = C' * R * C;

            % prefit residuals
            dY = Y_true(:, j) - Y_ACI;
            dy = C' * [dY(1:3);dY(5:6);dY(9)];

            % solve normal equations
            ax = H' * (r\H);
            nx = H' * (r\dy);
            
            Ax_L  = Ax_L + ax;
            Nx_L  = Nx_L + nx;
        end
    
        % solve LS
        XNOT_L = Ax_L\Nx_L;
    
        X_new(2:end) = X_new(2:end) + XNOT_L;
    
        [Cnm, Snm] = list2mat(n_max, Nc, Ns, X_new(1:Ncs));
    
        % update corrections
        xnot_L = xnot_L + XNOT_L;
    
        % show error
        err = vecnorm(XNOT_L);
        disp('  NSM update = '    + string(err));
    
        % update counter
        count = count + 1;
    end

    X_final = X_new;
    P_final = inv(Ax_L);
end

function [mxc, mcc] = compute_considerParams(R, Hc, Hn)
     mxc = (Hc' * (R \ Hn));
     mcc = (Hn' * (R \ Hn));
end

function [SK, SK2] = skewMatrix(vec)
    SK = -[0, vec(3), -vec(2);...
        -vec(3), 0 ,vec(1);...
         vec(2), -vec(1), 0];

    SK2 = SK*SK;
end