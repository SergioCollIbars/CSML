clear;
clc;
close all;

%% PLOT EART and MOON gravity field

% data text files
Earth_file  = 'GGM05S.txt';
Moon_file   = 'GRGM_400A.txt';

% Maximum gravity field
n_max  = 8;
Ncoeff = 0;
for j = 2:n_max
    Ncoeff = Ncoeff + (j + 1);
end

% load adata in matrix form
Earth_data = readmatrix(Earth_file);
Moon_data  = readmatrix(Moon_file);

% extract data 
C_earth  = Earth_data(1:Ncoeff, 4);
S_earth  = Earth_data(1:Ncoeff, 5);
Sg_c_earth = Earth_data(1:Ncoeff, 6);
Sg_s_earth = Earth_data(1:Ncoeff, 7);

C_moon  = Moon_data(3:2+Ncoeff, 3);
S_moon  = Moon_data(3:2+Ncoeff, 4);
Sg_c_moon = Moon_data(3:2+Ncoeff, 5);
Sg_s_moon = Moon_data(3:2+Ncoeff, 6);

% plot coefficients
lw = 2;
[num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
figure()
semilogy(1:Ncoeff, abs(C_earth), 'LineWidth', lw, 'color', 'r')
hold all;
semilogy(1:Ncoeff, abs(S_earth), 'LineWidth', lw, 'color', 'b')
semilogy(1:Ncoeff, abs(Sg_c_earth), 'LineWidth', lw, 'color', 'r', 'LineStyle', '--')
semilogy(1:Ncoeff, abs(Sg_s_earth), 'LineWidth', lw, 'color', 'b', 'LineStyle', '--')
title('Earth gravity field. Model GGM05C')
legend('C_{nm}','S_{nm}', '\sigma_C', '\sigma_S')
xticks(num_C);
xticklabels(str_C);

figure()
semilogy(1:Ncoeff, abs(C_moon), 'LineWidth', lw, 'color', 'r')
hold all;
semilogy(1:Ncoeff, abs(S_moon), 'LineWidth', lw, 'color', 'b')
semilogy(1:Ncoeff, abs(Sg_c_moon), 'LineWidth', lw, 'color', 'r', 'LineStyle', '--')
semilogy(1:Ncoeff, abs(Sg_s_moon), 'LineWidth', lw, 'color', 'b', 'LineStyle', '--')
title('Moon gravity field. Model GRGM660PRIM')
xticks(num_C);
xticklabels(str_C);

%% PLOT the tensor frobenius norm
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 
bodies = {'EARTH', '299', 'SUN', '5'}; % Earth, Venus, Sun , Jup Bar
et = cspice_str2et({'2020-01-10 12:00:00 TDB'});
ref = 'J2000';
abcorr = 'NONE';
observer = 'MOON';  % Set the observer to the Earth-Moon barycenter
for j = 1:length(bodies)
    % get GM
    [GM] = cspice_bodvrd(bodies{j}, 'GM', 1);                              % Get GM for the body [km^3/s^2]
    
    % set body
    target = bodies{j};
    [state, ~] = cspice_spkezr(target, et, ref, abcorr, observer);         % [Km & Km/s]
    posn = vecnorm(state(1:3));                                            % [Km]

    % compute tensor
    [T] = compute_tensor_GM(GM, posn);                                     % [1/s^2]

    % compute frobenius norm
    Tn = norm(T,"fro");                                                    % [1/s^2]
    Tn = Tn / (1E-9);                                                      % [Eotvos]

    disp("Body = " + bodies{j} + " Frobenius norm = " + string(Tn) + " E")
end

% clear kernels
cspice_kclear

%% FUNCTIONS
function [num_C, num_S, str_C, str_S] = SH_xlabel(n_max)    
    num_C = ones(1, n_max-1) * NaN;
    num_S = num_C;

    str_C = cell(1, n_max - 1);
    str_S = str_C;

    num_C(1) = 1;
    for j = 3:n_max
        num_C(j-1) = j + num_C(j-2);
    end
    
    num_S(1) = 1;
    for j = 3:n_max
        num_S(j-1) = (j-1) + num_S(j-2);
    end

    for j = 2:n_max
        str_C{j - 1} = "C_{" + string(j) + "0}";
        str_S{j - 1} = "S_{" + string(j) + string(j) + "}";
    end 
end

function [T] = compute_tensor_GM(GM, r)
    T = zeros(3,3);
    T(1, 1) = 2 * GM / r^3;     % [1/s^2]
    T(2, 2) = -GM / r^3;         % [1/s^2]
    T(3, 3) = -GM / r^3;         % [1/s^2]
end
