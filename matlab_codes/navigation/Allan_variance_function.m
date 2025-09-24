%% demo_allan_from_psd.m
% Example: white-noise signal -> PSD -> Allan deviation via integral

clear; clc;

%% 1) Make a test signal x(t)
Fs   = 100;           % Hz sampling
Tsec = 1800;          % 30 minutes
N    = Fs*Tsec;
t    = (0:N-1)'/Fs;

% White noise with std = 0.01 (units)
rng(1); 
x = 0.01 * randn(N,1);

%% 2) Estimate one-sided PSD Sx(f) with Welch
nfft   = 2^15;
wlen   = 8*Fs;                 % 8 s windows
nover  = round(0.5*wlen);      % 50% overlap
[Px, f] = pwelch(x, hamming(wlen), nover, nfft, Fs, 'onesided'); 
% Px units: (units^2 / Hz), f in Hz, f spans [0, Fs/2]

% If your PSD is of a DERIVATIVE and you want Allan of the BASE quantity,
% convert here. Example: given rate PSD S_omega(f), to get angle PSD S_theta(f):
% Px = Px ./ (2*pi*max(f,eps)).^2;   % uncomment if needed

%% 3) Choose averaging times and compute Allan deviation
tau = logspace(log10(1/Fs), log10(Tsec/10), 60);  % from ~1 sample to ~T/10
sigmaA = allan_from_psd(f, Px, tau);              % Allan deviation

%% 4) Plot
figure; 
loglog(tau, sigmaA, 'LineWidth', 2); grid on;
xlabel('\tau (s)');
ylabel('\sigma_A(\tau) (units)');
title('Allan Deviation from one-sided PSD');

% Optional: add a -1/2 slope guide (white noise reference)
hold on;
k = sigmaA(1)* (tau(1)^0.5);               % scale for visual guide
loglog(tau, k .* tau.^(-0.5), '--', 'LineWidth', 1);
legend('Allan deviation', 'slope -1/2 (white noise)', 'Location', 'southwest');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Helper function (same file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sigma = allan_from_psd(f, Sx, tau)
% f  : column vector of Hz (one-sided, starts at or > 0)
% Sx : same length as f, PSD of x(t) in (units^2/Hz)
% tau: vector of averaging times (s)
    f = f(:); 
    Sx = Sx(:);
    sigma2 = zeros(size(tau));
    for i = 1:numel(tau)
        x = pi.*f.*tau(i);
        W = (sin(x).^4) ./ (x.^2 + eps);   % Allan filter
        sigma2(i) = 2 * trapz(f, Sx .* W); % one-sided PSD -> factor 2
    end
    sigma = sqrt(sigma2);   % Allan deviation
end