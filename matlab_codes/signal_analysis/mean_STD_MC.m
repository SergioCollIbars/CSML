clear;
clc;
close all;
format long g;

%%                  MEAN, STD, MC SIMULATION
%   Description: Compute the mean and STD of a 1st order polyfit for the
%   1/f noise using the MC simulation


% Paths
addpath('functions/')

% Inputs
tmax = 86400;         % max simulation time [s]
N = 8640;             % Number of points time series
sigma_wn = 2.52E-10;   % STD GOCE (wn: White noise)
d = 6.32E-12 / sigma_wn;
sigma_pn = sigma_wn / (50 * d);   % PSD (pn: pink noise)

% Time and freq domains
x = linspace(1, tmax, N);
At = x(2) - x(1);
Fs = 1/At;
nfft = 2^nextpow2(N);
frec = 0 : Fs/(nfft-1) : Fs/2;

% define output values
p0 = ones(1, N) * NaN;
p1 = ones(1, N) * NaN;
STD = ones(1, N) * NaN;

% MC simulation
for j = 1:N
    disp("Iteration = " + string(j) + "/" +string(N));
    % generate noise
    [X, PSD] = noise_profile(sigma_wn, sigma_pn, N, Fs);
    
    % poyfit
    p = polyfit(x, X, 1);
    p0(j) = p(2);
    p1(j) = p(1);

    % Standart deviation
    STD(j) = std(diff(X)./At);
    clc
end

% compute resutls
p0_median = mean(p0);
p1_median = mean(p1);

p0_STD = std(p0);
p1_STD = std(p1);

STD_mean = mean(STD);

disp('Mean p0 = ' + string(p0_median))
disp('Mean p1 = ' + string(p1_median))

disp('STD p0 = ' + string(p0_STD))
disp('STD p1 = ' + string(p1_STD))

disp('STD time series = ' + string(STD_mean))

% plot histograms
figure()
histogram(p0)
title(" p0 histogram")

figure()
histogram(p1)
title("p1 histogram")

figure(STD)
title('STD histogram')



% % %% FUNCTION
% % function [Xhat] = LS(x, y)
% %     % reshape vectors
% %     N = length(x);
% %     x = reshape(x, [N, 1]);
% %     y = reshape(y, [N, 1]);
% % 
% %     H = [ones(N, 1), x];
% %     Ax = H' * H;
% %     Nx = H' * y;
% % 
% %     Xhat = Ax\Nx;
% % end