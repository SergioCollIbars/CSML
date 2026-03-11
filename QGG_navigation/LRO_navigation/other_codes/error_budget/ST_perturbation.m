clear; clc; close all;

set(0,'defaultAxesFontSize',16);
addpath('../../data/'); addpath(genpath("../../functions/"));
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels_LRO.tm')
%% SOLID TIDES PERTURBATION
% Description: Compute change in SH coefficients due to Solid Tides (ST).
% Date: 01/30/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GRGM1200 gravity field 
input_gravField  = "HARMCOEFS_MOON_GRGM1200.txt";

file             = readmatrix(input_gravField);
SH_coeff         = file(4:end);
[Nc, Ns, Ncs]    = count_num_coeff(file(1));

[Cnm, Snm] = list2mat(file(1), Nc, Ns, SH_coeff);

n2_static = [Cnm(3, 1:3),Snm(3,1:3)]';
n3_static = [Cnm(4, 1)];

%% Ephemerides (SPICE)
utc_start = '2012-03-01 00:00:00';
utc_stop  = '2012-05-01 00:00:00';
N         = 2000;                         % number of samples

et0 = cspice_str2et(utc_start);
et1 = cspice_str2et(utc_stop);
et  = linspace(et0, et1, N);

tgt       = 'SUN';
observer  = 'MOON';
ref_frame = 'MOON_PA'; % options: J2000 / IAU_MOON / MOON_PA
[sun_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

sun_SPICE(1:3, :) = sun_SPICE(1:3, :).*1E3;         % [m]
sun_SPICE(4:6, :) = sun_SPICE(4:6, :).*1E3;         % [m/s]

tgt       = 'EARTH';
observer  = 'MOON';
ref_frame = 'MOON_PA'; % options: J2000 / IAU_MOON / MOON_PA
[earth_SPICE, ~] = cspice_spkezr(tgt, et, ref_frame, 'NONE', observer);

earth_SPICE(1:3, :) = earth_SPICE(1:3, :).*1E3;         % [m]
earth_SPICE(4:6, :) = earth_SPICE(4:6, :).*1E3;         % [m/s]

% convert time to UTC
utc  = cspice_et2utc(et, 'ISOC', 6);
tUTC = datetime(utc, 'InputFormat', "yyyy-MM-dd'T'HH:mm:ss.SSSSSS");

%% Planet constants (SPICE)
[GM]       = cspice_bodvrd('MOON', 'GM', 1);  % Get GM for the Moon [km^3/s^2]
GM_M       = GM * 1E9; 
[radii]    = cspice_bodvrd('MOON', 'RADII', 3);
R_M        = radii(1).*1E3;                 % [m]

[GM]       = cspice_bodvrd('SUN', 'GM', 1);   % Get GM for the Sun [km^3/s^2]
GM_S       = GM * 1E9; 
[radii]    = cspice_bodvrd('SUN', 'RADII', 3);
R_S        = radii(1).*1E3;                  % [m]

[GM]       = cspice_bodvrd('EARTH', 'GM', 1); % Get GM for the Earth [km^3/s^2]
GM_E       = GM * 1E9; 
[radii]    = cspice_bodvrd('EARTH', 'RADII', 3);
R_E        = radii(1).*1E3;                   % [m]

%% Love Numbers from GRAIL
k20 = 0.02405;  s_k20 = 0.00018;
k21 = 0.02414;  s_k21 = 0.00025;
k22 = 0.02394;  s_k22 = 0.00028;
k30 = 0.0089;   s_k30 = 0.0021;

love_numb_list      = [k20, k21, k22, k30];
love_numb_list_unct = [s_k20, s_k21, s_k22, s_k30];

n_max               = 4;

%% Compute coefficient perturbation
MC      = 2;
n2_pert = nan(6, N, MC);
n3_pert = nan(1, N, MC);
for k = 1:N
    r_E = earth_SPICE(1:3, k);
    r_S = sun_SPICE(1:3, k);
    
    for mc = 1:MC
        [love_numb] = create_love_matrix(n_max, love_numb_list, ...
                    love_numb_list_unct);
        
        [deltaC, deltaS] = solid_tidal(n_max, n_max, love_numb, ...
            R_M, GM_M, GM_E, GM_S, r_E, r_S);

         n2_pert(:, k, mc) = [deltaC(3, 1:3),deltaS(3,1:3)]';
         n3_pert(:, k, mc) = [deltaC(4, 1)];
    end
end

% compute stats
N2_mean = mean(n2_pert, 3, 'omitnan');    % size: [6 x Nt]
N2_std  = std(n2_pert, 0, 3, 'omitnan');  % size: [6 x Nt]

N3_mean = mean(n3_pert, 3, 'omitnan');    % size: [6 x Nt]
N3_std  = std(n3_pert, 0, 3, 'omitnan');  % size: [6 x Nt]

figure; 
for i = 1:6
    mu = N2_mean(i,:);
    s3 = 3 * N2_std(i,:);

    % Mean line
    p = plot(tUTC, mu, 'LineWidth', 1.5); hold on;

    % Shaded region
    fill([tUTC', fliplr(tUTC')], ...
         [mu - s3, fliplr(mu + s3)], ...
         p.Color, ...
         'FaceAlpha', 0.4, ...
         'EdgeColor', 'none');
end
legend('\Delta C_{20}', '','\Delta C_{21}', '','\Delta C_{22}','', ...
    '\Delta S_{20}','', '\Delta S_{21}', '','\Delta S_{22}','');
grid on; title(' n = 2 perturbation')

figure; i = 1;
mu = N3_mean(i,:);
s3 = 3 * N3_std(i,:);

p = plot(tUTC, mu, 'LineWidth', 1.5); hold on;

fill([tUTC', fliplr(tUTC')], ...
     [mu - s3, fliplr(mu + s3)], ...
     p.Color, ...
     'FaceAlpha', 0.4, ...
     'EdgeColor', 'none');
legend('\Delta C_{30}'); title(' n = 3 perturbation'); grid on;


%% FUNCTIONS
function [deltaC, deltaS] = solid_tidal(n_max, n_sim, love_numb, ...
    R_M, GM_M, GM_E, GM_S, r_E, r_S)
    
    deltaC = zeros(n_max+1, n_max+1);  deltaS = zeros(n_max+1, n_max+1);

    % compute geographic coordinates
    phi_E    = atan2(r_E(3), (r_E(1)^2 + r_E(2)^2)^(0.5));
    phi_S    = atan2(r_S(3), (r_S(1)^2 + r_S(2)^2)^(0.5));
    lambda_E = atan2(r_E(2), r_E(1));
    lambda_S = atan2(r_S(2), r_S(1));
    r_E_norm = vecnorm(r_E);
    r_S_norm = vecnorm(r_S);

    % compute coefficient perturbation
    for n = 0:n_sim
        for m = 0:n
            row_index = n + 1; col_index = m+1;

            love_val  = love_numb(row_index, col_index);
            
            sinP_E = sin(phi_E);
            sinP_S = sin(phi_S);

            % Fully normalized associated Legendre
            Pbar   = legendre(n, sinP_E, 'norm');
            Pnm_E  = Pbar(m+1);

            Pbar   = legendre(n, sinP_S, 'norm');
            Pnm_S  = Pbar(m+1);

            Earth_pert = GM_E / GM_M * (R_M / r_E_norm)^(n+1) * Pnm_E;
            Sun_pert   = GM_S / GM_M * (R_M / r_S_norm)^(n+1) * Pnm_S;

            deltaC(row_index, col_index) = love_val / (2*n + 1) * ...
            (Earth_pert * cos(m * lambda_E) + Sun_pert * cos(m * lambda_S));

            deltaS(row_index, col_index) =  love_val / (2*n + 1) * ...
            (Earth_pert * sin(m * lambda_E) + Sun_pert * sin(m * lambda_S));
        end
    end
end

function [love_numb] = create_love_matrix(n_max, love_num_list, ...
    love_num_unct)
    love_numb = zeros(n_max+1, n_max+1);

    vals = love_num_list + love_num_unct .* randn(size(love_num_list));
     
    love_numb(3, 1:3)      = [vals(1), vals(2), vals(3)];
    love_numb(4, 1)        = vals(4);
end

