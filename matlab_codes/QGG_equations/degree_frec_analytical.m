clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);

%%                  FREQUENCY PER DEGREE ANALYTICAL
% Description: Compute the excited frequency for a certain degree in the
% gradiometer signal. 
%
% Date: 10/02/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference radius & planet mass
Ref = 1737.4E3;                     % [m] 
GM  = 4.902800118E12;               % [m^3/s^2] 
f_spin = 2*pi / (27.5 * 86400);     % [rad/s]

% load gravity coefficients
path = "HARMCOEFS_MOON_1_v2.txt";
[Cmat, Smat, ~] = readCoeff(path); % grav. field primary

% Mean orbital elements
e = 1E-2;              % [-]
i = deg2rad(90);    % [rad]
a = Ref + 90E3;     % [m]

% select mean orbit radius
r2_inv = 1/(a^2) * 1 / (sqrt(1-e^2));

% maximum degree
n_max = 1;

% compute mean secular rates (accounting only for J2 perturbation)
[Nnm] = nomalization_fnc(2, 0);
rates = j2_secular_rates(a, e, i, GM, Ref, -Nnm*Cmat(3,1));
Omega_dot = rates.Omega_rad_s.*0;
omega_dot = rates.omega_rad_s.*0;
M_dot     = rates.M_rad_s;

% create figure
figure; hold on;
xlabel('Frequency [Hz]');
ylabel('Amplitude [Eotvos]');
set(gca, 'XScale', 'log');
set(gca, 'YScale', 'log');
grid on;
title('Frequency-Amplitude for zonals');

colors = lines(n_max);
h_leg = []; % store one handle per color

for k = 0:n_max
    disp('Compute frequencies for degree ' + string(k) + ' ...');
    for m = 0:k
        % normalization for SH coefficients
        [Nnm] = nomalization_fnc(k, m);

        p_range = 0:1:k;
        q_range = 80;

        % go over p-values
        A = zeros(1, length(p_range)*q_range + 1);
        f = zeros(1, length(p_range)*q_range + 1);

        count = 1;
        for idx = 1:length(p_range)
            for q = -q_range/2:q_range/2
                p = p_range(idx);

                % compute inclination function
                F = compute_inclinationFunction(k, m, p, i);
    
                % compute frequency (2-side freq)
                f(count) = ((k - 2*p)*omega_dot + ...
                    (k - 2*p + q)*M_dot + m*(Omega_dot - f_spin))./(2*pi);
    
                % compute eccenticiy function
                G = compute_eccentricityFunction(k, p, q, e, 10);
                 
                % Amplitude per degree, order and p
                c = r2_inv * (k+1)*(k+2)*GM * ((Ref/ a)^k) * (1 / a);
                Cnm = Cmat(k+1, m+1); Snm = Smat(k+1, m+1);
                RMS_coeff = Nnm * sqrt(Cnm^2 + Snm^2);
    
                A(count) = 2*(c * (G * F) * RMS_coeff);
    
                % increment counter  
                count = count + 1;
            end
        end

        % select frequency range (only positive freq) & sum over common
        % frequencies
        [idx] = (f >= 0);
        f_pos = f(idx); A_pos = A(idx);
        [~, keep_idx, grp] = unique(f_pos, 'stable');% mapping to groups
        f_unique = f_pos(keep_idx);
        A_unique = A_pos(1:length(f_unique));
        
        for idx = length(f_unique)+1:length(f_pos)
            grpVal = grp(idx);
            val    = A_pos(idx);
            A_unique(grpVal) = A_unique(grpVal) + val; 
        end

        % Choose color
        c = colors(mod(k-1, size(colors,1)) + 1, :);

        % Draw vertical lines from 0 to amplitude
        fval = f_unique + 1E-12; Aval = abs(A_unique)./1E-9; % [Eotvos]
        ymin = 1e-10;
        if(m == 0)
            ls = '-';
            mk = 's';
        else
            ls = '-';
            mk = "diamond";
        end
        for d = 1:length(fval)
            plot([fval(d), fval(d)], [ymin, Aval(d)], ls, ...
                'Color', c, 'LineWidth', 1.5);
        end
    
        % Plot square markers at the top
         h = scatter(fval, Aval, 80, mk, 'filled', ...
            'MarkerFaceColor', c, 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
         if(m == 0)
            h_leg = [h_leg; h];  % append valid handle
         end
    end
end
legend_labels = arrayfun(@(n) sprintf('n = %d', n+1), 1:n_max-1, ...
    'UniformOutput', false);
legend(h_leg, legend_labels, 'Location', 'best');
ylim([1E-3, 10]);
xticks([1E-12 1E-5 1E-4 1E-3 1E-2, 1E-1])
xticklabels({'0','10^{-5}','10^{-4}','10^{-3}','10^{-2}', '10^{-1}'});

%% FUNCTIONS
function [Nnm] = nomalization_fnc(n, m)
 delta = 0;
 if(m == 0), delta = 1; end

 a = factorial(n-m); b = factorial(n+m);
 Nnm = sqrt((2-delta)*(2*n + 1) * a/b);
end
