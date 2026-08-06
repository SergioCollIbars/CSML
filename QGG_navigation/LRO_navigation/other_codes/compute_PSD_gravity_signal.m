clear;
clc;
close all;

set(0, 'defaultAxesFontSize', 16);
%% COMPUTE GRAVITY GRADIOMETER SIGNAL
% Author: Sergio Coll-Ibars
% Date: 07/09/2026


folder  = "/Users/sergiocollibars/Desktop/Lunar_orbit_simulator copy/";
file    = "LRO_traj_ideal.mat";

% load file
load(folder+file);

% reference frame
frame  = "RTN";
Y_new = Y;
if frame == "RTN"
    for j = 1:length(time)
        r = state(j, 1:3)'; v = state(j, 4:6)';
        [ECI_RTN] = RTN2ECI(r, v);
        RTN_ECI   = ECI_RTN';
        T_N       = [Y(1, j), Y(2, j), Y(3, j);...
                     Y(2, j), Y(4, j), Y(5, j);...
                     Y(3, j), Y(5, j), Y(6, j)];
        T_B       = RTN_ECI * T_N * RTN_ECI';
        Y_new(:, j) = [T_B(1,1);T_B(1,2);T_B(1,3);...
                       T_B(2,2);T_B(2,3);T_B(3,3)];

    end
end

% time increment [sec]
dt       = time(2) - time(1);
fs       = 1 /dt; % [Hz]
T_window = 3*60;

% time series
tt = ["XX", "XY", "XZ", "YY", "YZ", "ZZ"];
for k = 1:6
    figure();
    plot(time./3600, Y_new(k, :)/1E3, 'LineWidth', 1.5);
    title(tt(k)); grid on; xlabel('simulation time [hr]');
    ylabel('Eotvos');
end

% plot PSD signal
for k = 1:6
    figure();
    signal     = Y_new(k, :)/1E3;

    y_detrended = detrend(signal, 0); 

    N = length(y_detrended);

    segmentLength = min(4096, floor(N/4));

    segmentLength = max(segmentLength, 64);

    window     = hann(segmentLength);

    noverlap   = floor(0.5 * segmentLength);

    nfft       = max(4096, 2^nextpow2(segmentLength));
    [Pxx, f]   = compute_PSD(signal, fs, window, noverlap, nfft);
    [Pyy_full, f_full] = periodogram( ...
        y_detrended, hann(N), N, fs, 'onesided');

    ASD     = sqrt(Pxx);
    ASD_avg = abs(sinc(f*T_window)).*ASD;
    
    loglog(f_full(2:end), sqrt(Pyy_full(2:end)), ...
        'LineWidth', 0.5, 'Color', 'g'); hold on;
    loglog(f(2:end), ASD(2:end),     'LineWidth', 1.5,...
        'Color', 'b'); 
    loglog(f(2:end), ASD_avg(2:end), 'LineWidth', 1.5, ...
        'Color', 'k');
    
    title(tt(k)); grid on; xlabel('[Hz]'); ylabel('E Hz^{-0.5}');
    %ylim([1E-6 1E4]);
end

%% AUXILIARY FUNCTIONS
function [Pxx, f] = compute_PSD(X, Fs, window, noverlap, nfft)
    % plotPSD_multiChannel  Plot PSD for each channel (row) of X in one figure.
    %
    % INPUTS
    %   X        : N x Nt data matrix (N channels, Nt time samples)
    %   Fs       : sampling frequency [Hz]
    %   window   : (optional) pwelch window (scalar length or vector). Default: Hann(1024)
    %   noverlap : (optional) number of overlapping samples. Default: 50% of window length
    %   nfft     : (optional) FFT length. Default: max(256, 2^nextpow2(window length))
    %
    % NOTES
    %   - Uses Welch PSD estimate (pwelch)
    %   - Plots all channels on a single axes (log-log)
    %   - Assumes each row is a time series
    
    narginchk(2,5);
    
    if isempty(X) || ~ismatrix(X)
        error('X must be a non-empty 2D matrix of size N x Nt.');
    end
    if ~isscalar(Fs) || Fs <= 0
        error('Fs must be a positive scalar (Hz).');
    end
    
    [N, Nt] = size(X);
    
    % Defaults
    if nargin < 5 || isempty(window)
        winLen = min(1024, Nt);
        window = hann(winLen, 'periodic');
    else
        if isscalar(window)
            winLen = min(window, Nt);
            window = hann(winLen, 'periodic');
        else
            winLen = min(length(window), Nt);
            window = window(:);
            if length(window) ~= winLen
                window = window(1:winLen);
            end
        end
    end
    
    if nargin < 6 || isempty(noverlap)
        noverlap = floor(0.5 * winLen);
    else
        noverlap = min(noverlap, winLen-1);
    end
    
    if nargin < 7 || isempty(nfft)
        nfft = max(256, 2^nextpow2(winLen));
    end
    
    % Compute PSDs
    for i = 1:N
        xi = X(i,:);
        if all(~isfinite(xi))
            warning('Channel %d contains no finite values; skipping.', i);
            continue;
        end
    
        % Remove mean (helps with DC dominating)
        xi = xi - mean(xi, 'omitnan');
    
        % Replace NaNs/Infs (pwelch can't handle NaNs)
        bad = ~isfinite(xi);
        if any(bad)
            % simple linear interpolation over bad samples
            t = 1:Nt;
            good = ~bad;
            if nnz(good) < 2
                warning('Channel %d has too few valid samples; skipping.', i);
                continue;
            end
            xi(bad) = interp1(t(good), xi(good), t(bad), 'linear', 'extrap');
        end
    
        [Pxx, f] = pwelch(xi, window, noverlap, nfft, Fs, 'onesided');
    end
end

function [RTN_2_ECI] = RTN2ECI(r, v)
% ------------------------------------------------------------------- %
%                     RTN 2 ECI FRAME FUNCTION
% Author: Sergio Coll Ibars

% Date: 02/09/2023

% Description: Compute the radial, along-track, cross-track (RTN) 
%   to Earth-centered inertial rotation matrix. If applied to a
%   position vector in the RTN frame, it will transform that vector to
%   into the equivalent position vector in the ECI frame.

% Input:
%   r: current position vector. ECI frame
%   v: current velocity vector. ECI frame

% Output:
%   RTN_2_ECI: rotation matrix from RTN to ECI frame
% ------------------------------------------------------------------- %

n = cross(r, v);

R = r / vecnorm(r);
N = n / vecnorm(n);
T = cross(N, R);

RTN_2_ECI = [R, T, N];

end