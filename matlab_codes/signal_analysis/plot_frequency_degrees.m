clear;
clc;
close all;

%% COMPUTE FREQUENCY PER DEGREE
addpath("/Users/sergiocollibars/Desktop/CSML/QGG_gravEstim/data_files");
set(0,'defaultAxesFontSize',16);

orbit  = readtable('orbitData.txt');
time   = orbit.TIME;
state  = [orbit.ri_x,orbit.ri_y,orbit.ri_z]; % ACI frame [m]

path = "HARMCOEFS_BENNU_CD_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

n_array = 0:10;
for n_idx = 1:length(n_array)
    % current degree evaluated
    n       = n_array(n_idx);
    signal  = zeros(1, length(time)); 
    for k = 1:length(time)
        Wt = W * time(k) + W0;
        BN =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
        NB = BN';
        
        r_ECEF = BN * state(k, :)';
        
        % compute gravity. ECEF
        if(n > 0)
            [~, ~, ddU_ECEF_p] = potentialGradient_nm(Cnm, Snm, n, ...
                                                        r_ECEF, Re, GM, ...
                                                    normalized);
            [~, ~, ddU_ECEF_c] = potentialGradient_nm(Cnm, Snm, n-1, ...
                                                        r_ECEF, Re, GM, ...
                                                    normalized);
            ddU_ECEF = ddU_ECEF_c - ddU_ECEF_p;
        else
            [~, ~, ddU_ECEF] = potentialGradient_nm(Cnm, Snm, n, ...
                                                    r_ECEF, Re, GM, ...
                                                normalized);
        end
        ddU_ACI = NB * ddU_ECEF * BN;
        % signal(k)  = norm(ddU_ACI, 'fro');
        signal(k)  = ddU_ACI(3,3);
    end
    % compute periodogram
    fs = 1 / (time(2) - time(1));
    N  = length(time);
    X = fft(detrend(signal./1E-9));

    % two-sided
    P2 = (1/(N*fs)) * abs(X).^2;   % two-sided
    f2 = (0:N-1)*(fs/N);

    % one-sided
    if rem(N,2)==0
        P1 = P2(1:N/2+1);
        P1(2:end-1) = 2*P1(2:end-1);
        f1 = (0:N/2)*(fs/N);
    else
        P1 = P2(1:(N+1)/2);
        P1(2:end) = 2*P1(2:end);
        f1 = (0:(N-1)/2)*(fs/N);
    end
    % find peaks
    A = sqrt(P1);
    
    idx = f1 > 0 & isfinite(A) & A > 0;   % avoid zero frequency / log issues
    f = f1(idx);
    A = A(idx);

    if isempty(A)
        continue
    end
    % % [pks, locs] = findpeaks(A, f, 'MinPeakDistance', 5e-5);
    
    [pks, locs] = findpeaks(A, f, ...
    'MinPeakProminence', 0.0005*max(A), ...
    'MinPeakDistance', 1e-4);

    h = stem(locs, pks, 'filled');
    h.LineWidth = 1.5;
    h.MarkerSize = 4;

    % assign legend label
    h.DisplayName = "n = " + string(n);
    
    set(gca, 'XScale', 'log', 'YScale', 'log');
    hold on;
end
xlabel('Frequency [Hz]');
ylabel('$PSD [E / \sqrt{Hz}]$', 'Interpreter','latex');
grid on;
legend show;


