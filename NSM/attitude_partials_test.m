clear;
clc;
close all;

%% VERIFY ATTITUDE PARTIALS
% Verify the analytical atittude partials for the GG measurements.


% load planet parameters
[planetParams, poleParams, Kaula, ~, Xtrue, Iner] = loadPlanet("Earth");
GM  = planetParams(1); Re = planetParams(2); n_max = planetParams(3);
normalized = planetParams(4); [Nc, Ns, Ncs] = count_num_coeff(n_max); 

[Cnm, Snm] = list2mat(n_max, Nc, Ns, Xtrue);

% S/C location
r_hat = normrnd(0, 1, [3, 1]);
r_hat = r_hat./ vecnorm(r_hat);
r     = (Re + 250E3) * r_hat;  v = zeros(3, 1);

% compute true GG measurements
ACAF_ACI = eye(3); t = 0; state_t = [r', v'];
[Ytrue, ~, ~] = gradiometer_meas(t ,planetParams, ACAF_ACI, state_t, ...
                zeros(9, 1), Cnm, Snm);

T_ACI = [Ytrue(1), Ytrue(2), Ytrue(3);...
    Ytrue(4), Ytrue(5), Ytrue(6);...
    Ytrue(7), Ytrue(8), Ytrue(9)];

% body frame rotation
thy = pi/4; thp = pi/6; thr = pi/3; 
BN = rotationMatrix(thy, thp, thr, [3, 2, 1]);

% compute atittude error
yaw_err = .1E-6; pitch_err = .2E-6; roll_err = .3E-6;
RB = rotationMatrix(yaw_err, pitch_err, roll_err, [3, 2, 1]);

T_R = (RB*BN) * T_ACI * (RB*BN)';
T_B = BN * T_ACI * BN';

% compute difference
dT = reshape(T_R - T_B, [9, 1]);

% compute 1st order partials
[Hrot] = compute_rotPartials_analy(Ytrue, BN);
d      = Hrot * [yaw_err;pitch_err;roll_err];

% compute relative error
rel_err = abs(dT - d)./dT;