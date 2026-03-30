%% test_srif_recursive_demo.m
% Recursive SRIF demo: estimate constant parameters x from streaming linear measurements
% Plots estimation error and +/-3sigma bounds from SRIF covariance.

clear; clc; close all;
rng(7);  % reproducible

%% Problem setup
n  = 8;      % number of parameters to estimate
m  = 6;      % measurements per epoch
Nt = 200;    % number of time steps

% True parameter vector
x_true = randn(n,1);

% Measurement noise (example: diagonal covariance, can be full SPD too)
sigma = 0.05;                     % measurement noise std (same for all m)
Rk = (sigma^2) * eye(m);          % measurement covariance

% Storage
xhat_hist = nan(n,Nt);
Pdiag_hist = nan(n,Nt);
err_hist  = nan(n,Nt);

%% SRIF initialize (no a priori)
Rinfo = zeros(n,n);
zinfo = zeros(n,1);

%% Main loop (recursive SRIF)
for k = 1:Nt
    % Make a "reasonably exciting" Hk so the problem becomes observable
    % (Random + slowly varying structure)
    Hk = randn(m,n) + 0.2*sin(2*pi*k/30)*randn(m,n);

    % Simulate measurement
    vk = chol(Rk,'lower')*randn(m,1);
    yk = Hk*x_true + vk;

    % Residual form: y = H x + v  -> treat residual vector as y (since model is linear)
    % SRIF wants equation H*x ≈ y in whitened LS sense.
    rk = yk;

    % --- SRIF measurement update for this epoch
    [Rinfo, zinfo] = srif_meas_update(Rinfo, zinfo, Hk, rk, Rk);

    % Solve estimate and covariance at each step
    x_hat = Rinfo \ zinfo;            % backsolve (R is upper triangular)
    P_hat = srif_cov_from_R(Rinfo);   % P = inv(R'*R)

    % Store
    xhat_hist(:,k) = x_hat;
    Pdiag_hist(:,k) = diag(P_hat);
    err_hist(:,k) = x_hat - x_true;
end

%% Plot: error and +/-3sigma bounds for a few parameters
idx = 1:min(4,n);    % plot first 4 parameters (change as you like)
t = 1:Nt;

figure;
for i = 1:numel(idx)
    j = idx(i);
    subplot(numel(idx),1,i);
    sigj = sqrt(Pdiag_hist(j,:));

    plot(t, err_hist(j,:), 'LineWidth', 1.5); hold on;
    plot(t,  3*sigj, '--', 'LineWidth', 1.2);
    plot(t, -3*sigj, '--', 'LineWidth', 1.2);
    grid on;
    ylabel(sprintf('e_%d', j));
    if i==1
        title('SRIF Recursive Estimation Error with \pm 3\sigma Bounds');
    end
    if i==numel(idx)
        xlabel('time step k');
    end
    legend('error', '+3\sigma', '-3\sigma', 'Location', 'best');
end

%% Plot: RMS error norm vs time and a "typical" 3-sigma norm bound
err_norm = sqrt(sum(err_hist.^2,1));
sig_norm = 3*sqrt(sum(Pdiag_hist,1));   % crude bound using diag(P) only

figure;
plot(t, err_norm, 'LineWidth', 1.8); hold on;
plot(t, sig_norm, '--', 'LineWidth', 1.3);
grid on;
xlabel('time step k');
ylabel('||e||_2');
title('Parameter Error Norm vs 3\sigma Diagonal Bound');
legend('||x̂-x||', '3\sigma \approx 3\sqrt{trace(P)}', 'Location', 'best');

disp('Done. If the error stays mostly within ±3σ, your SRIF + covariance convention is consistent.');

%% -------------------- Local functions --------------------

function [Rinfo_new, zinfo_new] = srif_meas_update(Rinfo, zinfo, Hk, rk, Rk)
% One recursive SRIF measurement update:
% Given prior (Rinfo, zinfo) representing Rinfo*x ≈ zinfo,
% incorporate measurement equation Hk*x ≈ rk with covariance Rk.

    n = size(Rinfo,1);

    % Whitening using Cholesky (lower): Rk = S*S'
    S = chol(Rk,'lower');

    Htil = S \ Hk;   % = inv(S)*Hk  (whitened)
    rtil = S \ rk;   % = inv(S)*rk

    % Stack and QR
    A = [Rinfo, zinfo;
         Htil,  rtil];

    [~, Raug] = qr(A,0);

    Rinfo_new = Raug(1:n, 1:n);
    zinfo_new = Raug(1:n, end);

    % Enforce positive diagonal (optional, makes covariance consistent)
    d = diag(Rinfo_new);
    sgn = ones(size(d));
    sgn(d<0) = -1;
    D = diag(sgn);

    Rinfo_new = D*Rinfo_new;
    zinfo_new = D*zinfo_new;
end

function P = srif_cov_from_R(Rinfo)
% Compute P = inv(R'*R) stably, without forming inverse explicitly.
    n = size(Rinfo,1);
    Ri = Rinfo \ eye(n);
    P = Ri*Ri';
end
