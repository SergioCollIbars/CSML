clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);

format long g;
%%                 FISHER INFORMATION BASED ON SPECTRAL LINES
% Description: Compute the excited frequency for a certain degree in the
% gradiometer signal. Then, compute the Signal to Noise Ratio SNR
%
% Date: 10/13/2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reference radius & planet mass (Moon)
Ref = 1737.4E3;                     % [m] 
GM  = 4.902800118E12;               % [m^3/s^2] 
f_spin = 2*pi / (27.5 * 86400);     % [rad/s]

% load gravity coefficients
path = "HARMCOEFS_MOON_1_v2.txt";
[Cmat, Smat, ~] = readCoeff(path); % grav. field primary

n_max = 20;

% Mean orbital elements
i = 0:1:90;
i = deg2rad(i);          % [rad]
e = 0:1E-2:0.1;          % [-]
a = Ref + 90E3;          % [m]
T = 2*pi / sqrt(GM/ a^3); % [sec]

% Information matrix per inclination
SNR = ones(length(i), length(e)) * NaN;

% Measurement band
MB = [1E-4, 1E-1];

for incIdx = 1:length(i)
    disp('Computing inclination value: ' + string(rad2deg(i(incIdx))))
    for eccIdx = 1:length(e)
        % select inclination
        inc = i(incIdx);
        ecc = e(eccIdx);
    
        % compute amplitude & frequency
        [Amplitude_nm,frequency_nm] = compute_Amplitude_Freq(n_max, GM, Ref,...
            ecc, inc, a, Cmat, Smat, f_spin);
        
        % compute total fisher information matrix
         [F, ~, ~] = SNR_from_lines(frequency_nm, Amplitude_nm, ...
           T, MB, 1E-10);
    
        SNR(incIdx, eccIdx)  = F;
    end
end

% plot results
[x, y] = meshgrid(e, rad2deg(i));
figure()
contourf(x, y, 10*log10(SNR))
colorbar()
xlabel('Eccentricity [-]')
ylabel('Inclination [deg]')
title('SNR in dB for \Gamma_{rr}')
% % grid on;
% % xlabel('Inclination [deg]');
% % ylabel('Information');

%% FUNCTIONS
function [Nnm] = nomalization_fnc(n, m)
 delta = 0;
 if(m == 0), delta = 1; end

 a = factorial(n-m); b = factorial(n+m);
 Nnm = sqrt((2-delta)*(2*n + 1) * a/b);
end

function [Amplitude_nm,frequency_nm] = compute_Amplitude_Freq(n_max, ...
    GM, Ref, e, i, a, Cmat, Smat, f_spin)

    r2_inv = 1/(a^2) * 1 / (sqrt(1-e^2));
    
    number_coeff = 0;
    for k = 2:n_max
        number_coeff = number_coeff + (k + 1);
    end

    % compute mean secular rates (accounting only for J2 perturbation)
    [Nnm] = nomalization_fnc(2, 0);
    rates = j2_secular_rates(a, e, i, GM, Ref, -Nnm*Cmat(3,1));
    Omega_dot = rates.Omega_rad_s;
    omega_dot = rates.omega_rad_s;
    M_dot     = rates.M_rad_s;

    % go over p-values
    Amplitude_nm = ones(number_coeff , 1E4)*NaN;
    frequency_nm = ones(number_coeff , 1E4)*NaN;
    count2 = 1;
    for k = 2:n_max
        for m = 0:k
            % normalization for SH coefficients
            [Nnm] = nomalization_fnc(k, m);
    
            p_range = 0:1:k;
            q_range = 4;
    
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
                    c = r2_inv * (k+1)*(k+2)*(k+3)*GM * ((Ref/ a)^k) * (1 / a);
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
    
            % save amplitudes & frequencies per degree order
            Amplitude_nm(count2, 1:length(A_unique)) = A_unique;
            frequency_nm(count2, 1:length(f_unique)) = f_unique;
    
            count2 = count2 + 1;
        end
    end
    
    % Eliminate NaNs
    [mask] = ~isnan(frequency_nm);
    valid_cols = all(mask, 1);
    frequency_nm = frequency_nm(:, valid_cols);
    Amplitude_nm = Amplitude_nm(:, valid_cols);
end


function [F, f_bins, A_sum_per_bin] = SNR_from_lines(f_mat, A_mat, T_orb, band, df)
% f_mat, A_mat : size [nLines x nFreqCols] (same shape). Each row is (n,m); cols are frequencies.
% T_orb        : observation time (one orbit) [s]
% band         : [fmin fmax] Hz (e.g., [1e-4 1e-2])
% df           : bin width Hz (use df = 1/T_orb if using DFT binning)
% Spsd         : function handle for noise PSD, Eotvos^2/Hz, e.g., @(f) S0
%
% Returns:
%   F              : SNR for the single scalar parameter
%   f_bins         : center frequency per bin (Hz)
%   A_sum_per_bin  : sum of amplitudes in each bin (Eotvos/m)

    arguments
        f_mat double
        A_mat double
        T_orb (1,1) double {mustBePositive}
        band (1,2) double = [0, inf]
        df (1,1) double {mustBePositive} = 1/max(T_orb,eps)
    end

    % 1) Flatten and mask invalids
    f = f_mat(:);
    A = A_mat(:);
    ok = ~isnan(f) & ~isnan(A);
    f = f(ok);  A = A(ok);

    % 2) Band-pass selection
    fmin = min(band); fmax = max(band);
    inband = (f >= fmin) & (f <= fmax);
    f = f(inband);  A = A(inband);

    if isempty(f)
        F = 0; f_bins = []; A_sum_per_bin = [];
        return
    end

    % 3) Bin frequencies (tolerance via df). Use nearest-bin indexing.
    %    This is robust to tiny floating point differences.
    bin_idx = round(f / df);   % integers
    % Sum amplitudes per bin:
    A_sum_per_bin = accumarray(bin_idx, A, [], @sum);

    % specify measurement sigma / sqrt(Hz)
    sigma = 1E-12;

    % 4) Recover bin-center frequencies and keep only bins that appeared
    used_bins = find(A_sum_per_bin ~= 0);
    A_sum_per_bin = A_sum_per_bin(used_bins);
    f_bins = used_bins * df;

    % 5) Fisher information (single parameter): F = (T/2) * sum ( (sum A)^2 / S(f) )
    F = sum((A_sum_per_bin.^2)) / (sigma^2);
end


