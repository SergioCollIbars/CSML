clear;
clc;
close all;
addpath('functions/');
set(0,'defaultAxesFontSize',16);

%%                  ATTITUDE ERROR PARAMETRIZATION
% Description: Evaluate the attitude residual for different magnitudes of
% attitude errors by using the NSM in Earth applications.
% Author: Sergio Coll-Ibars
% Date: 04-25-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                Mesurement mask 
%           xx xy xz yx yy yz zx zy zz
mask     =  [1, 1, 1, 0, 1, 1, 0, 0, 1]';
dir      =  [2 3]; Ndir = length(dir);
names    = ["xx", "xy", "xz", "yx", "yy", "yz", "zx", "zy", "zz"];
%% Extract GG observations
folderPath = "/Users/sergiocollibars/Documents/GG_observations/120by120/500km/";
[GG_obs]   = parser_GG_obs_MSODP(folderPath); Nt = length(GG_obs(:, 1));

%% Attitude error magnitudes (arcseconds)
Att_noise_norm = logspace(-3,3, 100); Natt = length(Att_noise_norm);
scale          = 3600 * 180 / pi; % [arcseconds/rad]

%% Compute residuals
att_res_norm_NSM = nan(Ndir, Natt);
att_res_norm_ORG = nan(sum(mask), Natt);
for k = 1:Natt
    disp(k/Natt.*100);
    att_Err = normrnd(0, Att_noise_norm(k)/scale, [3, Nt]);
    [GG_rot, GG_nom]  = rotate_GG_obs(GG_obs,att_Err);

    prefit            = GG_rot(:, 2:end) - GG_nom(:, 2:end);
    dY                = prefit(:, logical(mask))'; % mask x Nt
    dY_NSM            = nan(Ndir, Nt);
    for j = 1:Nt
        [Hrot_1st]    = compute_rotPartials_analy(GG_nom(j, 2:end)',...
            eye(3));
        [S1, V1, D1]  = svd(Hrot_1st(logical(mask), :)');
        V_rot         = D1(:, dir);
        dY_NSM(:, j)  = V_rot' * dY(:, j);
    end
    for g = 1:Ndir
        att_res_norm_NSM(g, k)   = rms(dY_NSM(g, :))./1E-9; % [Eotvos]
    end
    for p = 1:sum(mask)
        att_res_norm_ORG(p, k)   = rms(dY(p, :))./1E-9;     % [Eotvos]
    end
end

%% Plot results
figure()
loglog(Att_noise_norm, att_res_norm_NSM, 'LineWidth', 1.2, 'Color', 'k');
grid on; xlabel('arcseconds'); ylabel('Eotvos');
title('RMS NSM residuals');

figure()
loglog(Att_noise_norm, att_res_norm_ORG, 'LineWidth', 1.5);
grid on; xlabel('arcseconds'); ylabel('Eotvos');
title('RMS Original residuals');
legend(names(logical(mask)));