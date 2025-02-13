clear;
clc;
close all;

%%      MEASUREMENT EFFECT IN GM ERROR

% ADD PATH
addpath('functions/')

% TIME INPUT
N = 1000;
tmin = 0;
tmax = 86400 * 5;
TIME = linspace(tmin, tmax, N);

% SH matrices
Cnm_t = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.0391557863539988 -2.96928723209235e-06 0.00375640654748657 0 0 0 0;...
        0.0148427177700986 0.00167095097673949 3.80845003468165e-05 0.000371755938456641 0 0 0;...
        0.0307494000000000 0.000413625917950024 -0.000490123739988179 -6.43092753252371e-05 4.51227856599555e-05 0 0;...
        -0.00456599734888228 0.000441961635589938 8.84264903209225e-05 1.40396087936725e-05 1.07458402471302e-05 1.19886311472876e-06 0;...
        -0.00896736657720649 0.000905916675449317 9.05780703119059e-05 2.77025499573633e-05 -1.92680576794137e-06 9.97533032070693e-08 1.67034838314692e-07];

Snm_t = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 0 -2.54325906400954e-05 0 0 0 0;...
        0 0.000992000134408593 6.53000000000000e-05 -0.00100797120329237 0 0 0;
        0 0.000634013000392474 0.000108054642426876 0.000102400000000000 0.00291093983173820 0 0;...
        0 -1.75754943031483e-05 -7.41878397813858e-05 4.79413750994751e-06 -0.000503800000000000 0.000448812426298560 0;...
        0 -1.35941163743731e-06 -9.69889209840526e-06 -7.55736000569855e-06 -1.59676058718751e-06 -0.000599300000000000 -3.93397896234722e-05];

% ORBIT INPUTS
GM = 5.2;           % Gravitational parameter [m^3/s^2]
Re = 246;           % Planet reference radius [m]
a = 2000;           % semi-major axis [m]
e = 0.5;           % eccectricity [-]
i = deg2rad(0);     % inclination [rad]
omega = 0;          % argument of periapsis [rad]
Omega = 0;          % RAAN [rad]
f = 0;              % true anomaly [rad]
rho = a*(1 - e*e);  % orbital parameter [m]
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% inital state. IJK frame
r0 = rho / (1 + e * cos(f)) * [cos(f); sin(f); 0];
v0 = sqrt(GM / rho) * [-sin(f); e + cos(f); 0];
[BN] = rotationMatrix(Omega, i, omega, [3,1,3]);
r0 = BN' * r0;
v0 = BN' * v0;
X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];


% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, GM, 0, Re, 1), TIME, X0, options);

r_ACI = state(:, 1:3)';
v_ACI = state(:, 4:6)';


% MEASUREMENTS
deltaGM = 2;
errRTN = zeros(N, 9);
for j = 1:N

     Wt = W0 + W * TIME(j);
     ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

    r = r_ACI(:, j);
    v = v_ACI(:, j);
    [ACI_RTN] = RTN2ECI(r, v);


    [~, ~, ddU1] = potentialGradient_nm(Cnm_t, Snm_t, 1, ...
                                           ACAF_ACI*r, Re, GM, 1);
    ddU1 = ACAF_ACI' * ddU1 * ACAF_ACI;

    [~, ~, ddU2] = potentialGradient_nm(Cnm_t, Snm_t, 1, ...
                                           ACAF_ACI*r, Re, GM+deltaGM, 1);
     ddU2 = ACAF_ACI' * ddU2 * ACAF_ACI;

    err = ACI_RTN' * abs(ddU1 - ddU2) * ACI_RTN;
    errRTN(j, :) = reshape(err, [9, 1]);
end 

