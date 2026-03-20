clear;
clc;
close all;
addpath('../../NSM/functions/EoM_functions/')
addpath('../../QGG_gravEstim/')
addpath('../../QGG_gravEstim/src/')
addpath('../../QGG_gravEstim/data_files/')
addpath('../../QGG_navigation/data/')
set(0,'defaultAxesFontSize',16);

%%                  DEGREE VARIANCE EQUATION
% Description: plot the theorical gravity estimation accuracy for the
% gradiometer-based gravity estimation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Asteroid parameters.
path = "HARMCOEFS_BENNU_CD_1.txt";
[Cnm, Snm, R] = readCoeff(path);
GM = 5.2;
n_max  = 16;
normalized = 1;
n_range = [0, 2:n_max];

% orbit radius
% r = 1E3;                 % [m]
r_norm = 4; 
r  = r_norm * R;         % [m]
disp('radius: ' + string(r));
disp('Altitude: ' + string(r - R));
T = sqrt(r^3 / GM);      % [sec] 

% Extract gravity information
[Nc, Ns, Ncs] = count_num_coeff(n_max); 
[X] = mat2list(Cnm, Snm, Nc, Ns);
X_RMS = computeRMS_coeffErr(n_max, Nc, Ns, ...
            X, zeros(n_max+1, n_max+1), zeros(n_max+1, n_max+1));

asterParams = [GM, R, n_max, normalized];

% sampling frequency
f = 1/10;              % [Hz]
t_max = 9 * 86400;     % [sec]
L = t_max * f;         % [-]
rev = t_max / T;
disp('Number of meas: ' + string(L));

% measurement accuracy 
sigma_m = 1E-12;       % [s^-2]

RMS_theorical = ones(1, length(n_range));
for k = 1:length(n_range)
    n = n_range(k);
    lambda = (n+1) * (n+2);
    
    RMS_theorical(k) = r^3/GM * (r/R)^n * sigma_m / sqrt(L) * ...
        1 / sqrt(lambda);
end

% plot 
figure()
semilogy(n_range(2:end), RMS_theorical(2:end), 'LineWidth', 1.5, ...
    'Marker','*', 'Color', 'b', 'LineStyle','-')
hold on;
semilogy(n_range(2:end), X_RMS(2:end), 'Color', 'k', 'LineWidth', 2)
xlabel('Degree [n]')
ylabel('RMS [-]')
grid on;