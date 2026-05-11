clear; clc;
close all;
format long g;
addpath('functions/');

%% COMPUTE ANGULAR VELOCITY from the nominal attitude obs.
% Files generated with MSODP
% Author: Sergio Coll-Ibars
% Date: 03/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FOLDER = "/Users/sergiocollibars/Documents/att_files/410km";
files  = dir(FOLDER);
nfiles = length(files);

% save angular velocity, options: 1 or 0
save_angVel = 1;

quaternions = nan(1e8, 4); idx = 1;
for k = 1:length(files)
    file = files(k).name;
    if(contains(file, ".att"))
        disp(files(k).name)

        att_noise_free  = read_att_to_rpy(FOLDER+'/'+file);

        sod = att_noise_free.sod;
        if(sod(1)~=86385 || sod(17284)~=0)
            a = 1;
        end
        vals = 4:17283;
        data = att_noise_free.q(vals, :);

        N = length(data);
        quaternions(idx:N+idx-1, :) = data;
        idx = idx+N;
    end
end

% time delta
dt = att_noise_free.t(2) - att_noise_free.t(1); % [s]

% eliminate NANs
[idx] = ~isnan(quaternions(:, 1));
q       = quaternions(idx, :);

% time derivative
q_dot_0 = gradient(q(:, 1), dt);
q_dot_1 = gradient(q(:, 2), dt);
q_dot_2 = gradient(q(:, 3), dt);
q_dot_3 = gradient(q(:, 4), dt);
q_dot   = [q_dot_0, q_dot_1, q_dot_2, q_dot_3];

% compute angular velocity
angVel  = zeros(3, length(q(:, 1)));
for k = 1:length(q(:, 1))
    q0 = q(k, 1);
    q1 = q(k, 2);
    q2 = q(k, 3);
    q3 = q(k, 4);

    E = [-q1, q0, q3, -q2;...
         -q2, -q3, q0, q1;...
         -q3, q2, -q1, q0];

    angVel(:, k) = 2 * E * q_dot(k, :)';
end

% remove outliers
omega1 = filloutliers(angVel(1, :), "pchip", "movmedian", 101);
omega2 = filloutliers(angVel(2, :), "pchip", "movmedian", 101);
omega3 = filloutliers(angVel(3, :), "pchip", "movmedian", 21);

angVel = [omega1;omega2;omega3];
angAcc = nan(3, length(angVel(1, :)));
for k  = 2:length(angVel(1, :))-1
    angAcc(:, k) = (angVel(:, k+1) - angVel(:, k-1))./(2*dt);
end

figure();

Nt = 1:1:length(q(:, 1));
t = Nt.*5 - 5;
plot(t, angVel, 'LineWidth', 2); grid on;
title('Angular velocity'); ylabel('[rad/s]');
legend('\omega_1', '\omega_2', '\omega_3');

figure();

plot(t, angAcc, 'LineWidth', 2); grid on;
title('Angular acceleration'); ylabel('[rad/s^2]');
legend('\Omega_1', '\Omega_2', '\Omega_3');

if (save_angVel == 1)
    disp('saving angular velocity');
    save("nominal_angVel_410km.mat", 'angVel');
end