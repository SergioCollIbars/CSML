clear; clc; close all;

set(0,'defaultAxesFontSize',16);
addpath('../../data/'); addpath(genpath("../../functions/"));
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm');

%% STUDY STATE CORRELATION
% Description: Study correlation between the states
% Date: 01/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GRGM1200 gravity field 
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);
[Nc, Ns, Ncs]    = count_num_coeff(file(1));

[Cnm, Snm] = list2mat(file(1), Nc, Ns, SH_coeff);

normalized = file(3);

%% mask
%       xx xy xz yy yz zz
mask = [1, 1, 1, 1, 1, 1];


%% Ephemerides (SPICE)
utc_start = '2025-03-16 00:00:00';
utc_stop  = '2025-03-17 00:00:00';
N         = 200;                              % number of samples
[GM]      = cspice_bodvrd('MOON', 'GM', 1);   % Get GM for the Moon [km^3/s^2]
GM_M      = GM * 1E9;                         % [m^3/s^2]
[radii]   = cspice_bodvrd('MOON', 'RADII', 3);
R_M       = radii(1).*1E3;                    % [m]

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

tgt       = 'LUNAR RECONNAISSANCE ORBITER';
observer  = 'MOON';
ref_frame = 'J2000'; % options: J2000 / IAU_MOON
[sc_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

sc_SPICE(1:3, :) = sc_SPICE(1:3, :).*1E3;         % [m]
sc_SPICE(4:6, :) = sc_SPICE(4:6, :).*1E3;         % [m/s]

%% Compute SC and planet orientation
attitude  = "RTN"; % options: inertial / RTN

BN_MOON_mat = nan(3 * N, 3); BN_mat = nan(3 * N, 3);
for k = 1:N
    maxInd = 3 * k; minInd = maxInd - 2;
    
    % S/C orientation
    if(attitude == "inertial")
        BN = eye(3);
    elseif(attitude == "RTN")
        NB = RTN2ECI(sc_SPICE(1:3, k), sc_SPICE(4:6, k));
        BN = NB';
    end
    BN_mat(minInd:maxInd, :) = BN;

    % planet orientation
    frame_to   = 'J2000';
    frame_from = 'MOON_PA';
    J2000_MOON = cspice_pxform(frame_from, frame_to, et(k));
    BN_MOON_mat(minInd:maxInd, :)  = J2000_MOON';

end
%% Correlations
traj      = ones(6, 1); att = ones(6, 1).*0; bias = ones(6, 1); 
stateMask = [traj;att;bias];

n_max     = 150;
SF        = ones(6, 1); 

sigma_GG  = 1;                       % [milli-Eotvos]
sigma_att = 1 * pi / (180 * 3600);   % [rad]
W_GG      = eye(6) * ((1/sigma_GG)^2);
W_GGM     = eye(sum(mask)) * ((1/sigma_GG)^2);
W_att     = eye(3) * ((1/sigma_att)^2);

W         = diag([diag(W_GGM); diag(W_att)]);
% % W         = W_GG; 

Ax_min     = zeros(6,6); Ax_max = Ax_min;
states_std = nan(6, N);
figure(); labels = {'R','T','N','\simga_1','\sigma_2', '\sigma_3'}; % example names
h =  imagesc(zeros(6));          % Plot matrix as image
colorbar;            % Show color scale
axis equal tight;    % Keep square aspect
colormap(jet);       % Choose colormap (optional)

xticks(1:length(labels));
yticks(1:length(labels));

xticklabels(labels);
yticklabels(labels);
for k = 1:N
    maxInd = 3 * k; minInd = maxInd - 2;
    BODYMOON_J2000 = BN_MOON_mat(minInd:maxInd, :);

    r = sc_SPICE(1:3, k); v = sc_SPICE(4:6, k); % J2000 frame
    [Y_J2000, ~] = gradiometer_meas(et(k) ,[GM_M, R_M, n_max, normalized],...
        BODYMOON_J2000, [r', v'], zeros(9, 1), Cnm, Snm);

    % rotation to Instrument frame
    BN = BN_mat(minInd:maxInd, :);

    T_ACI = [Y_J2000(1),Y_J2000(2),Y_J2000(3);...
             Y_J2000(4),Y_J2000(5),Y_J2000(6);...
             Y_J2000(7),Y_J2000(8),Y_J2000(9)];
    T_B   = BN * T_ACI * BN';

    Y     = [T_B(1,1);T_B(1,2);T_B(1,3);T_B(2,2);...
               T_B(2,3);T_B(3,3)]./1E-12;   % [mE]

    % compute measurement partials
    BODYMOON_BODY = BODYMOON_J2000 * BN'; r_ACI = r;
    [Hpos] = compute_posPartials(n_max, normalized, Cnm, Snm, R_M, GM_M, ...
                r_ACI, BODYMOON_J2000, BODYMOON_BODY);
    Hp     = SF.*[Hpos(1:3, :); Hpos(5:6, :);Hpos(9, :)]./1E-12;

    % compute attitude partials
    [Hrot] = compute_rotPartials_analy(Y_J2000, BN);
    Hr     = SF.*[Hrot(1:3, :); Hrot(5:6, :);Hrot(9, :)]./1E-12;

    G      = Hp' * W_GG * Hr;

    np = sqrt(diag(Hp.' * W_GG * Hp));    % 3x1, weighted norms of Hr columns
    nt = sqrt(diag(Hr.' * W_GG * Hr));    % 3x1, weighted norms of Hth columns

    cosA = G ./ (np * nt.');              % 3x3 matrix of cos(alpha_ij)

    maxCorr = max(max(abs(cosA)));
    minCorr = min(min(abs(cosA)));

% %     H       = [Hp, Hr];
    H       = [Hp(logical(mask), :), Hr(logical(mask), :);...
               zeros(3), eye(3)];

    Ax_min = Ax_max; 
    Ax_max = Ax_min + H' * W * H; 

    if(rank(Ax_max) == 6)
        P = inv(Ax_max);

        states_std(:, k)   = diag(sqrt(P));
        corr_mat = corr(P);
        cosA_filt = corr_mat(1:3, 4:6);
        disp(max(max(abs(cosA_filt))));
    else
        corr_mat = zeros(6);
    end
    % Update plot data
    set(h, 'CData', abs(corr_mat));
    drawnow;  % Force update
end


figure()
semilogy(et, vecnorm(states_std(1:3, :)), 'LineWidth', 2);
title('Position std'); grid on;

figure()
semilogy(et, vecnorm(states_std(4:6, :)), 'LineWidth', 2);
title('Attitude std'); grid on;
