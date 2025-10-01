
function ekf_augmented_bias_demo()
% EKF with state-augmented random-walk bias (1D constant-velocity example)
% - Truth generation (position, velocity, bias)
% - Additive measurement with bias
% - Analytic and numerical measurement partials (Jacobian)
% - EKF estimation with 3-sigma plots & consistency checks

%% -------------------- Settings --------------------
rng(7);                 % Reproducible
dt   = 30;             % [s] sample time
T    = 2000;             % [s] total time
Nt   = round(T/dt);     % steps

% Continuous-time white accel PSD driving kinematics (pos/vel model)
Sa   = 1e-2;            % [m^2/s^3] (tune for roughness of motion)

% Continuous-time random-walk PSD for bias
Sb   = 1e-1;            % [unit^2/s] (tune for bias drift rate)

% Measurement noise covariance
R    = (0.5)^2;         % [unit^2] measurement variance

% Initial truth and filter guesses
x0_true   = [0; 1.0];   % [pos; vel]
b0_true   = 2.5;        % bias initial value (units of measurement)
x0_hat    = [0; 0.0];   % filter initial state guess
b0_hat    = 0.0;        % filter initial bias guess

P0_x   = diag([10^2, 1^2]); % initial cov on [pos; vel]
P0_b   = (1.0)^2;           % initial cov on bias

%% -------------------- Models --------------------
% State: x = [pos; vel]; bias b is random-walk.
% Continuous model: xdot = [0 1; 0 0] x + [0; 1] * a_w,  a_w ~ N(0, Sa)
F  = [0 1; 0 0];
G  = [0; 1];

% Discrete kinematics (exact for CV model)
Phi_x = [1 dt; 0 1];
Qx    = Sa * [ (dt^3)/3, (dt^2)/2; (dt^2)/2, dt ];    % standard CV integral

% Bias: b_{k+1} = b_k + w_b, w_b ~ N(0, Sb*dt)
Phi_b = 1;
Qb    = Sb*dt;

% Augmented discrete transition and process noise
Phi = blkdiag(Phi_x, Phi_b);
Q   = blkdiag(Qx, Qb);

% Measurement: z = h(x,b) + v = pos + b + v
% Analytic measurement Jacobian (partials):
% H = [d h / d pos, d h / d vel, d h / d b] = [1, 0, 1]
h_fun  = @(x_aug) x_aug(1) + x_aug(3);
H_fun  = @(x_aug) [1, 0, 1];  %#ok<NASGU> (kept for clarity)

% Optional numerical Jacobian checker
jac_check = false;  % set true to verify analytic H against numeric
if jac_check
    H_num = numjacobian(h_fun, [x0_hat; b0_hat]);
    disp("Numeric H at initial guess:");
    disp(H_num);
end

%% -------------------- Generate Truth & Measurements --------------------
x_true = zeros(2, Nt);
b_true = zeros(1, Nt);
z      = zeros(1, Nt);
t      = (0:Nt-1)*dt;

x_true(:,1) = x0_true;
b_true(1)   = b0_true;

for k = 1:Nt-1
    % Process noise draws
    w_x = mvnrnd([0;0], Qx).';
    %w_b = sqrt(Qb) * randn;
    w_b = normrnd(0, sqrt(Qb), [1,1]);

    % Propagate truth
    x_true(:,k+1) = Phi_x * x_true(:,k) + w_x;
    b_true(k+1)   = Phi_b * b_true(k)   + w_b;
end

% Measurements
v = sqrt(R) * randn(1, Nt);
for k = 1:Nt
    z(k) = x_true(1,k) + b_true(k) + v(k);
end

%% -------------------- EKF Initialization --------------------
xhat = zeros(3, Nt);
Phat = zeros(3, 3, Nt);

xhat(:,1) = [x0_hat; b0_hat];
Phat(:,:,1) = blkdiag(P0_x, P0_b);

% Pre-allocate storage for stats
innov   = zeros(1, Nt);
S_store = zeros(1, Nt);
NEES    = zeros(1, Nt);
NIS     = zeros(1, Nt);

%% -------------------- EKF Loop --------------------
H_analytic = [1 0 1];

for k = 2:Nt
    % ---------- Predict ----------
    xhat_pred = Phi * xhat(:,k-1);
    P_pred    = Phi * Phat(:,:,k-1) * Phi.' + Q;

    % ---------- Update ----------
    % (Using analytic H; you can swap for numeric if needed)
    H = H_analytic;
    h = xhat_pred(1) + xhat_pred(3);

    y = z(k) - h;                       % innovation
    S = H * P_pred * H.' + R;           % innovation cov
    K = (P_pred * H.') / S;             % Kalman gain

    xhat(:,k)   = xhat_pred + K*y;      % posterior
    Phat(:,:,k) = (eye(3) - K*H) * P_pred * (eye(3) - K*H).' + K*R*K.'; % Joseph

    % ---------- Consistency metrics ----------
    innov(k)   = y;
    S_store(k) = S;

    err = [x_true(:,k); b_true(k)] - xhat(:,k);
    NEES(k) = err.' / Phat(:,:,k) * err;   % normalized estimation error squared
    NIS(k)  = y / S * y;                   % normalized innovation squared
end

%% -------------------- Plots --------------------
figure('Name','State & Bias Estimation','Color','w'); 
tiledlayout(3,1, 'Padding','compact', 'TileSpacing','compact');

% Position
nexttile;
plot(t, x_true(1,:), 'k-', 'LineWidth', 1.5); hold on;
plot(t, x_true(1,:) - xhat(1,:), 'LineWidth', 1.2);
sig_pos = squeeze(sqrt(Phat(1,1,:))).';
plot(t, +3*sig_pos, '--', t, -3*sig_pos, '--');
ylabel('pos [m]'); grid on; legend('true','error','\pm3\sigma','Location','best');

% Velocity
nexttile;
plot(t, x_true(2,:), 'k-', 'LineWidth', 1.5); hold on;
plot(t, x_true(2,:) - xhat(2,:), 'LineWidth', 1.2);
sig_vel = squeeze(sqrt(Phat(2,2,:))).';
plot(t, +3*sig_vel, '--', t, -3*sig_vel, '--');
ylabel('vel [m/s]'); grid on; legend('true','error','\pm3\sigma','Location','best');

% Bias
nexttile;
plot(t, b_true, 'k-', 'LineWidth', 1.5); hold on;
plot(t, b_true - xhat(3,:), 'LineWidth', 1.2);
sig_b = squeeze(sqrt(Phat(3,3,:))).';
plot(t, +3*sig_b, '--', t, -3*sig_b, '--');
ylabel('bias [unit]'); xlabel('time [s]'); grid on;
legend('true','error','\pm3\sigma','Location','best');

figure('Name','Innovations & Consistency','Color','w');
tiledlayout(3,1, 'Padding','compact', 'TileSpacing','compact');

% Innovations
nexttile;
plot(t, innov, 'LineWidth', 1); hold on;
plot(t, 3*sqrt(S_store), 'r--', t, -3*sqrt(S_store), 'r--');
ylabel('innovation'); grid on; legend('y','\pm3\sqrt{S}');

% NIS
nexttile;
plot(t, NIS, 'LineWidth', 1);
ylabel('NIS'); grid on;

% NEES
nexttile;
plot(t, NEES, 'LineWidth', 1);
ylabel('NEES'); xlabel('time [s]'); grid on;

disp('Done. Tip: tune Sa, Sb, and R to explore consistency and bias tracking.');
end

%% ---------- Utility: numerical Jacobian (finite differences) ----------
function J = numjacobian(hfun, x, eps)
% NUMJACOBIAN  Finite-difference Jacobian for scalar h(x) at x
if nargin < 3 || isempty(eps), eps = 1e-6; end
n = numel(x);
J = zeros(1,n);
h0 = hfun(x);
for i = 1:n
    dx = zeros(n,1); dx(i) = eps;
    J(i) = (hfun(x+dx) - h0)/eps;
end
end
