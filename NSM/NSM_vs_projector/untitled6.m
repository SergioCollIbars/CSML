clear;
clc;
close all;

% Asteroid / planet parameters
GM      = 5.2;   % [Kg m^3/s^2]
n_max   = 0;

path    = "HARMCOEFS_BENNU_CD_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
asterParams = [GM, Re, n_max, 1];
% S/C orbit
y = Re + 600;
z = 0;
x = 0;

% rotation states
rn_ACI      = [x;y;z];
ACAF_ACI    = eye(3);
ACI_BODY    = eye(3);
ACAF_BODY   = ACAF_ACI * ACI_BODY;

% compute position partials
[Hpos] = compute_posPartials(n_max, 1, Cnm, Snm, Re, GM, ...
                rn_ACI, ACAF_ACI, ACAF_BODY);

Hn = [Hpos(1:3, :);Hpos(5:6, :);Hpos(9,:)];
disp(Hn);
disp(rank(Hn));

% compute attitude partials
[Y_ACI, ~] = gradiometer_meas(0 ,asterParams, ACAF_ACI,...
                [rn_ACI', zeros(3, 1)'], ...
                    zeros(9,1), Cnm, Snm);
[Hatt] = compute_rotPartials_analy(Y_ACI, ACI_BODY');
 Hn = [Hatt(1:3, :); Hatt(5:6, :);Hatt(9, :)];
disp(Hn);
disp(rank(Hn));