clear;
clc;
close all;
format long g;
addpath('../../NSM/functions/')
addpath('../../QGG_gravEstim/src/')
addpath('GOCE_L1b_MatlabReaders_1.1/data/')
addpath('GOCE_L2b_MatlabReaders/data/')
set(0,'defaultAxesFontSize',16);

% Earth planet data
path = "HARMCOEFS_EARTH_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
path = "SIGMACOEFS_EARTH_1.txt";
[sigma_Cnm, sigma_Snm, ~] = readCoeff(path);
GM = 3.986004418E14;
n_max  = 20;
normalized = 1;
W = 0;                      % Rotation ang. vel   [rad/s]
W0 = 0;                     % Initial asteroid longitude
RA = -pi/2;                 % Right Ascension     [rad]
DEC = pi/2;                 % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);

saveData = 1;

% trajectory data
date = "16_Nov_2012";
% % var = ["pos","vel","posCov","velCov", "time"];
% % for j = 1:5
% % file = var(j)+"_"+date+".mat";
% % load(file);
% % end
load(date + "_L2position.mat");
[commonTimes, idx1, idx2] = intersect(TT_GPS_PVT_final, positions(:, 1));
pos_L1 = POS_PVT_FINAL(idx1, :);
pos_L2 = positions(idx2, 2:end);
Nt     = length(commonTimes);
Nt     = round(Nt/1); t = commonTimes(1:Nt);
pos_L2 = pos_L2(1:Nt, :);
pos_L1 = pos_L1(1:Nt, :);

Ar = 1.5;                                           % position rms error [m]

% generate gradiometer data
sigma = 1E-12;                                      % [1/s^2]
noise = normrnd(0, sigma, [9, Nt]); noise0 = zeros(9, 1);
R = diag([sigma, sigma, sigma, sigma, sigma].^2);   % meas. weight
if(saveData)
    filename = "gradMeas.mat";
    if ~isfile(filename)  % or use exist(filename, 'file') == 2
        [Ytrue, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, pos_L2, ...
                noise, Cnm, Snm, eye(3,3));
        save("gradMeas.mat", "Ytrue");
    else
        disp('File exists. Loading...');
        load(filename);
    end
else
    [Ytrue, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, pos_L2, ...
                noise, Cnm, Snm, eye(3,3));
end

% run estimation
[S] = mat2list(sigma_Cnm, sigma_Snm, Nc, Ns);
sigma_RMS     = computeRMS_coeffErr(n_max, Nc, Ns, S, Cnm.*0, Snm.*0);
sigma_n = 100*sigma_RMS(2:end);
[Xp, Pp] = perturb_coeff(sigma_n, n_max, X);
[Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);
P0 = Pp(2:end, 2:end); 

% loop
iterMax = 1;
count   = 0;
xnot_N = zeros(Ncs-1, 1); xnot_LS = xnot_N;
Cp_N = Cp; Cp_LS = Cp;
Sp_N = Sp; Sp_LS = Sp;
Xp_N = Xp; Xp_LS = Xp;

 % apriori values for the Consider Parameters; 
Pc = zeros(3, 3); Pxc = zeros(Ncs - 1, 3); Pxc_N = zeros(Ncs - 1, 6);
Pc_N = zeros(6, 6);
for j = 1:length(Pc)
    Pc(j, j) = (Ar)^2;
end
for j = 1:length(Pc_N)
    Pc_N(j, j) = Ar(1)^4;
end
c = zeros(3, 1).*Ar;
c_N = zeros(6, 1).*Ar(1)^2;

while count < iterMax
    Ax_N = inv(P0); Ax = inv(P0);
    Nx_N = -inv(P0) * xnot_N; Nx = -inv(P0) * xnot_LS;

    [~, Mxc, Mcc] = get_considerCov_apriori(P0, Pc, Pxc);
    [~, Mxc_N, Mcc_N] = get_considerCov_apriori(P0, Pc_N, Pxc_N);
    h = waitbar(0, 'Processing...');
    for j = 1:Nt
         waitbar(j / Nt, h, sprintf('Processing: %d of %d', j, Nt));
        % current position
        rn_ECEF   = pos_L1(j, :);
        rn_ECEF2  = pos_L2(j, :);
        ACAF_ACI  = eye(3,3);

        if(prod(~isnan(Ytrue(:, j))) && prod(~isnan(rn_ECEF)))
            % Null space merthod + Apriori (AP)
            [Yc, Hc_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn_ECEF, zeros(1, 3)], ...
                    noise0, Cp_N, Sp_N, eye(3,3));
    
            [Hpos] = compute_posPartials(n_max, normalized, Cp_N, Sp_N, Re, GM, rn_ECEF', ACAF_ACI, ACAF_ACI);
            
            Hap = compute_posPartials_2ndOrder(GM, rn_ECEF(1), rn_ECEF(2), rn_ECEF(3));
            [ax, nx, mxc, mcc] = nullSpace_method_AP(Ytrue(:, j)-Yc, Hc_ACI, R, Hpos, Hap, eye(3,3), zeros(9, 1));
            Ax_N  = Ax_N  + ax;
            Nx_N  = Nx_N  + nx;
            Mxc_N = Mxc_N + mxc;
            Mcc_N = Mcc_N + mcc;
            
            % LS method
            [Yc, Hc_ACI, ~] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn_ECEF, zeros(1, 3)], ...
                    noise0, Cp_LS, Sp_LS, eye(3,3));
    
            [Hpos] = compute_posPartials(n_max, normalized, Cp_LS, Sp_LS, Re, GM, rn_ECEF', ACAF_ACI, ACAF_ACI);
    
            [ax, nx, mxc, mcc] = LS_method(Ytrue(:, j)-Yc, Hc_ACI, R, Hpos, eye(3,3), zeros(9, 1));
            Ax = Ax + ax;
            Nx = Nx + nx;
            Mxc = Mxc + mxc;
            Mcc = Mcc + mcc;
        end
    end
    close(h)

    % solve LS
    XNOT_LS = Ax\Nx - Ax\(Mxc * c);
    XNOT_N = Ax_N\Nx_N - Ax_N\(Mxc_N * c_N);

    Xp_N(2:end) = Xp_N(2:end) + XNOT_N;
    Xp_LS(2:end) = Xp_LS(2:end) + XNOT_LS;

    [Cp_N, Sp_N] = list2mat(n_max, Nc, Ns, Xp_N);
    [Cp_LS, Sp_LS] = list2mat(n_max, Nc, Ns, Xp_LS);

    % update corrections
    xnot_N = xnot_N + XNOT_N;
    xnot_LS = xnot_LS + XNOT_LS;

    % show error
    disp('Null space update = ' + string(vecnorm(XNOT_N)));
    disp('LS update = ' + string(vecnorm(XNOT_LS)));
    disp('--------------------------------------------------------');           

    % update counter
    count = count + 1;
end
Px = inv(Ax);
Sxc = -Px * Mxc;
Pxx = Px + Sxc*Pc*Sxc';
Pxc = Sxc * Pc;

Px_N = inv(Ax_N);
Sxc_N = -Px_N * Mxc_N;
Pxx_N = Px_N + Sxc_N*Pc_N*Sxc_N';
Pxc_N = Sxc_N * Pc_N;

P_LS =  [Pxx, Pxc;Pxc', Pc];
P_N =  [Pxx_N, Pxc_N;Pxc_N', Pc_N];
sigma_LS = sqrt(diag(P_LS));
sigma_N = sqrt(diag(P_N));

[Xp_N] = mat2list(Cp_N, Sp_N, Nc, Ns);
[Xp_LS] = mat2list(Cp_LS, Sp_LS, Nc, Ns);

SH_N = Xp_N(2:end);
SH_LS = Xp_LS(2:end);

% compute RMS error & plot
X_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        X, Cnm.*0, Snm.*0); 
val = [1; sqrt(diag(P0))];
sigma0_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
       val, Cnm.*0, Snm.*0);  
err_init_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        X - Xp, Cnm.*0, Snm.*0);
errN_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        X - [0;SH_N], Cnm.*0, Snm.*0);
errLS_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
        X - [0;SH_LS], Cnm.*0, Snm.*0);
figure()
semilogy(1:n_max, X_RMS, 'LineWidth', 2, 'Color', 'k')
grid on;
hold all;
semilogy(1:n_max, sigma0_RMS, 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
 semilogy(1:n_max, err_init_RMS, 'LineWidth', 2, 'Color', 'r')
semilogy(1:n_max, errN_RMS, 'LineWidth', 2)
semilogy(1:n_max, errLS_RMS, 'LineWidth', 2)
legend('truth', '\sigma 0', 'init error', 'NSM level 1', 'LS level 2')
title('RMS order error')
xlabel('order')
xlim([2, n_max])

%% FUNCTIONS
function [ax, nx, mxc, mcc] = nullSpace_method_AP(Y, Hc, R, Hp, Hap, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;

    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];

    % look for null space
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');

    % project measurements
    y  = C' * dY;
    hc = C' * [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
        Hc(8, 2:end)];
    hap = C' * Hap;
    r  = C' * R * C;

    % information and normal matrices
    ax = hc' * inv(r) * hc;
    nx = hc' * inv(r) * y;
    mxc = (hc' * inv(r) * hap);
    mcc = (hap' * inv(r) * hap);
end

function [ax, nx, mxc, mcc] = LS_method(Y, Hc, R, Hp, ACI_RTN, noise)
        % reshape meas [ACI]
        ddU = [Y(1), Y(2), Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
        % Rotate meas to RTN coordinates
        dy_RTN = ACI_RTN' * ddU * ACI_RTN;
        dy = reshape(dy_RTN, [9, 1]) + noise;

        hc = [Hc(1, 2:end); Hc(4, 2:end); Hc(7, 2:end);Hc(5, 2:end);...
                Hc(8, 2:end)];
       
        hp = [Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)];
    
        % select measurements
        dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];
 
        ax  = hc' * inv(R) * hc;
        nx  = hc' * inv(R) * (dY);
        mxc = (hc' * inv(R) * hp);
        mcc = (hp' * inv(R) * hp); 
end
