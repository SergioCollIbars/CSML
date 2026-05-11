clear; clc;
close all;
format long g;
addpath('functions/');

%% COMPUTE ATTITUDE NOISE from the nominal vs true observation attitude 
% Files generated with MSODP
% Author: Sergio Coll-Ibars
% Date: 03/13/2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FILE_noise_free = "/Users/sergiocollibars/Documents/att_files/" + ...
    "grace.2008-08-01_A_61_nonoise.att.pre";
att_noise_free  = read_att_to_rpy(FILE_noise_free);
DCM_noise_free  = read_att_to_DCM(FILE_noise_free);

FILE_noise = "/Users/sergiocollibars/Documents/att_files/" + ...
    "grace.2008-08-01_A_61_mode4.att.pre";
att_noise = read_att_to_rpy(FILE_noise);
DCM_noise = read_att_to_DCM(FILE_noise);

Nt      = length(DCM_noise(1, 1, :));
DCM_err = nan(3,3, Nt);
for k = 1:Nt
    DCM_err(:,:, k) =  DCM_noise_free(:,:,k)' * DCM_noise(:, :, k); % convetion in MSODP: R_actual = R_nom * R_noise
end

err_rpy = nan(3, Nt); % yaw-pitch-roll (3-2-1)
for k =1:Nt
    C = DCM_err(:,:, k);
    y = atan2(C(1,2), C(1,1)); % yaw
    p = -asin(C(3,1));         % pitch
    r = atan2(C(2,3), C(3,3)); % roll

    err_rpy(:, k) = [r;p;y];
end

% compute error PSD
dt = 5; % seconds
[f, PSD_r] = compute_PSD(err_rpy(1, :), dt);
[~, PSD_p] = compute_PSD(err_rpy(2, :), dt);
[~, PSD_y] = compute_PSD(err_rpy(3, :), dt);

% compute angular velocity (nominal)
dt = 5;
angVel = nan(3, Nt);
for k = 2:Nt-1
    DCM_prev = DCM_noise_free(:, :, k-1);
    DCM_curr = DCM_noise_free(:, :, k);
    DCM_next = DCM_noise_free(:, :, k+1);

    DCM_dot = (DCM_next - DCM_prev) ./ (2 * dt);

    omega_dyad = - DCM_dot *(DCM_curr');
     
    angVel(1, k) = (omega_dyad(3, 2) - omega_dyad(2, 3))/2;
    angVel(2, k) = (omega_dyad(1, 3) - omega_dyad(3, 1))/2;
    angVel(3, k) = (omega_dyad(2, 1) - omega_dyad(1, 2))/2;
end

figure()
plot(1:Nt, angVel); grid on; title('Angular velovity body frame');
% % save('nominal_angVel.mat', 'angVel');
%% plot time series. Noise
figure;
plot(att_noise.t, err_rpy(1, :))
xlabel('Time since first sample [s]')
ylabel('Roll [rad]')
grid on

figure;
plot(att_noise.t, err_rpy(2, :))
xlabel('Time since first sample [s]')
ylabel('Pitch [rad]')
grid on

figure;
plot(att_noise.t, err_rpy(3, :))
xlabel('Time since first sample [s]')
ylabel('Yaw [rad]')
grid on

%% plot PSD Noise
figure(); lw = 1.2;
loglog(f, PSD_r, 'LineWidth', lw); hold on;
loglog(f, PSD_p, 'LineWidth', lw); 
loglog(f, PSD_y, 'LineWidth', lw); 
legend('roll', 'pitch', 'yaw'); grid on;
title('Attitude error PSD'); xlabel('Hz'); ylabel('[rad^2 / Hz]')

%% plot yaw-pitch-roll orientation w.r.t inertial frame
figure;
plot(att_noise.t,  att_noise_free.roll)
xlabel('Time since first sample [s]')
ylabel('Roll [rad]')
grid on

figure;
plot(att_noise.t,att_noise_free.pitch)
xlabel('Time since first sample [s]')
ylabel('Pitch [rad]')
grid on

figure;
plot(att_noise.t, att_noise_free.yaw)
xlabel('Time since first sample [s]')
ylabel('Yaw [rad]')
grid on

RMS_roll  = rms(err_rpy(1, :));
RMS_pitch = rms(err_rpy(2, :));
RMS_yaw   = rms(err_rpy(3, :));

norm_att_err = vecnorm(err_rpy);
RMS_att_err  = rms(norm_att_err);
MAX_att_err  = max(norm_att_err);

arcseconds  = 180 * 3600 / pi; % [arc/ rad]
disp('RMS attitude error: ');
disp(RMS_att_err * arcseconds);

disp('MAX attitude error: ');
disp(MAX_att_err * arcseconds);

%% PLOT TIME FRAME EVOLUTION
% roll, pitch, yaw are Nx1 vectors
% time is Nx1
% angles assumed in radians

N = Nt;
time  = att_noise_free.t;

% Initial state @ 07/31/2008
P0= [1758567.79573676     111794.11238752    6640146.72986794];
V0= [-7360.62324564081      50.97471707645    1952.23388529411];

x0 = [P0, V0];

muE = 398600.4418E9; % km^3/s^2

opts = odeset('RelTol',1e-12,'AbsTol',1e-12);

[~, state] = ode113(@(t,x) two_body_eom(t,x,muE), time, x0, opts);
xout = state(:, 1:3);

figure; hold on; grid on; axis equal;
xlabel('X'); ylabel('Y'); zlabel('Z');
view(3);

% Trajectory
plot3(xout(:,1), xout(:,2), xout(:,3), 'k--', 'LineWidth', 1.0);

% Scale arrow length relative to trajectory size
traj_scale = max(vecnorm(xout - mean(xout,1), 2, 2));
L = 0.05 * traj_scale;

% Inertial frame arrows at origin
quiver3(0,0,0,L,0,0,'k','LineWidth',2);
quiver3(0,0,0,0,L,0,'k','LineWidth',2);
quiver3(0,0,0,0,0,L,'k','LineWidth',2);

text(L,0,0,'I_X');
text(0,L,0,'I_Y');
text(0,0,L,'I_Z');

% Initial spacecraft position
r0_plot = xout(1,:);

% Initial body-frame arrows located at spacecraft position
hBx = quiver3(r0_plot(1), r0_plot(2), r0_plot(3), L,0,0,'r','LineWidth',2);
hBy = quiver3(r0_plot(1), r0_plot(2), r0_plot(3), 0,L,0,'g','LineWidth',2);
hBz = quiver3(r0_plot(1), r0_plot(2), r0_plot(3), 0,0,L,'b','LineWidth',2);

% Optional marker for spacecraft position
hSc = plot3(r0_plot(1), r0_plot(2), r0_plot(3), 'ko', ...
            'MarkerFaceColor','k', 'MarkerSize', 6);

% Axis limits based on trajectory
pad = 0.1 * traj_scale;
xlim([min(xout(:,1))-pad, max(xout(:,1))+pad]);
ylim([min(xout(:,2))-pad, max(xout(:,2))+pad]);
zlim([min(xout(:,3))-pad, max(xout(:,3))+pad]);

roll  = att_noise_free.roll;
pitch = att_noise_free.pitch;
yaw   = att_noise_free.yaw;

N = min([length(roll), length(pitch), length(yaw), length(time), size(xout,1), 1e3]);

for k = 1:N

    % Current spacecraft position
    r_sc = xout(k,:);

    phi   = roll(k);
    theta = pitch(k);
    psi   = yaw(k);

    R1 = [1 0 0;
          0 cos(phi) -sin(phi);
          0 sin(phi)  cos(phi)];

    R2 = [ cos(theta) 0 sin(theta);
           0          1 0;
          -sin(theta) 0 cos(theta)];

    R3 = [cos(psi) -sin(psi) 0;
          sin(psi)  cos(psi) 0;
          0         0        1];

    % If this is body-to-inertial:
    R = R3 * R2 * R1;

    % Body-frame axes expressed in inertial coordinates
    bx = R(:,1);
    by = R(:,2);
    bz = R(:,3);

    % Update arrow origins and directions
    set(hBx, ...
        'XData', r_sc(1), 'YData', r_sc(2), 'ZData', r_sc(3), ...
        'UData', L*bx(1), 'VData', L*bx(2), 'WData', L*bz(3)*0 + L*bx(3));

    set(hBy, ...
        'XData', r_sc(1), 'YData', r_sc(2), 'ZData', r_sc(3), ...
        'UData', L*by(1), 'VData', L*by(2), 'WData', L*by(3));

    set(hBz, ...
        'XData', r_sc(1), 'YData', r_sc(2), 'ZData', r_sc(3), ...
        'UData', L*bz(1), 'VData', L*bz(2), 'WData', L*bz(3));

    % Update spacecraft marker
    set(hSc, ...
        'XData', r_sc(1), ...
        'YData', r_sc(2), ...
        'ZData', r_sc(3));

    title(sprintf('Time = %.2f s', time(k)));

    drawnow;
end
%% FUNCTIONS
function dxdt = two_body_eom(~, x, mu)

    r = x(1:3);
    v = x(4:6);

    a = -mu * r / norm(r)^3;

    dxdt = [v; a];

end