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
thy = 0; thp = 0; thr = 0; 
BN = rotationMatrix(thy, thp, thr, [3, 2, 1]);

% compute atittude error
yaw_err = .1E-8; pitch_err = .2E-8; roll_err = .3E-8;
RB = rotationMatrix(yaw_err, pitch_err, roll_err, [3, 2, 1]);

T_R = (RB*BN) * T_ACI * (RB*BN)';
T_B = BN * T_ACI * BN';

% compute difference
dT = reshape(T_R - T_B, [9, 1]);

% compute 1st order partials
[Hrot] = compute_rotPartials_analy(Ytrue, BN);
[H] = rotation_partials_v2(T_ACI);
d      = Hrot * [yaw_err;pitch_err;roll_err];

% compute relative error
rel_err = abs(dT - d)./dT;



%% FUNCTIONS
function [H] = rotation_partials_v2(T)
    S1 = [0,0,0;0,0,1;0,-1,0];
    S2 = [0,0,-1;0,0,0;1,0,0];
    S3 = [0,1,0;-1,0,0;0,0,0];

    dT_d1 = S1 * T + (S1 * T)';
    dT_d2 = S2 * T + (S2 * T)';
    dT_d3 = S3 * T + (S3 * T)';

    h_d1 = [dT_d1(1,1);dT_d1(1,2);dT_d1(1,3);dT_d1(2,1);...
        dT_d1(2,2);dT_d1(2,3);dT_d1(3,1);dT_d1(3,2);dT_d1(3,3)];
    h_d2 = [dT_d2(1,1);dT_d2(1,2);dT_d2(1,3);dT_d2(2,1);...
        dT_d2(2,2);dT_d2(2,3);dT_d2(3,1);dT_d2(3,2);dT_d2(3,3)];
    h_d3 = [dT_d3(1,1);dT_d3(1,2);dT_d3(1,3);dT_d3(2,1);...
        dT_d3(2,2);dT_d3(2,3);dT_d3(3,1);dT_d3(3,2);dT_d3(3,3)];

    H = [h_d3, h_d2, h_d1];
   
end