clear;
clc;
close all;
format long g;
set(0,'defaultAxesFontSize',16);

%%              PLOT ORBITAL ELEMENTS
% Description: Plot the orbital elements for the selected conditions
% Author: Sergio Coll
% Date: 10/27/2025

target_body = "Eros";  % options: Bennu or Eros

if(target_body == "Bennu")
    % Asteroid parameters.
    path = "HARMCOEFS_BENNU_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    GM = 5.2;
    n_max  = 6;
    normalized = 1;
    W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
    W0 = 0;                   % Initial asteroid longitude
    RA = deg2rad(86.6388);    % Right Ascension     [rad]
    DEC = deg2rad(-65.1086);  % Declination         [rad]
    r = 0.36E3;               % orbit radius        [m]
elseif(target_body == "Eros")
    path = "HARMCOEFS_EROS_CD_1.txt";
    [Cnm, Snm, Re] = readCoeff(path);
    n_max  = 4;
    normalized = 1;
    GM =  459604.431484721;          % Point mass value    [m^3/s^2]
    W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
    W0 = 0;                          % Initial asteroid longitude
    RA = deg2rad(11.363);            % Right Ascension     [rad]
    DEC = deg2rad(17.232);           % Declination         [rad]
    r = 24E3;                        % orbit radius        [m]
end
poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n   = sqrt(GM / r^3);    % Mean motion         [rad/s]
T   = (2 * pi / n);
rev = 3;
f   = 1/1;
t   = linspace(0, rev*T, rev*T * f);
Nt  = length(t);

% Integrate trajectory (true trajectory)
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, 14, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0;STM0], options);


% compute orbital elements
orbital = ones(6, Nt) * NaN;
for k = 1:Nt
    r_ACI = state_t(k, 1:3)'; v_ACI = state_t(k, 4:6)'; 
    [alpha]= orbitalElem(r_ACI, v_ACI, asterParams(1));
    orbital(:, k) = [alpha(1);alpha(3);rad2deg(wrapToPi(alpha(5))); ...
        rad2deg(alpha(6));rad2deg(alpha(7));...
        rad2deg(wrapToPi(alpha(8)))];
end

figure()
tt = ["e [-]", "a [m]", "f [deg]", "i [deg]", "RAAN [deg]", ...
    "\omega [deg]"];
for k = 1:length(orbital(:, 1))
    subplot(2, 3, k)
    plot(t./T, orbital(k, :), 'LineWidth', 2, 'Color', 'b')
    grid on;
    xlabel('Orbital rev.')
    ylabel(tt(k))
    if(k == 3 || k == 6)
        ylim([-180 180]);
        yticks([-180 -100 0 100 180]);
    end
end
sgtitle('Orbital Elements ' + string(target_body))

init_orbital  = orbital(:, 1);
disp('Initial orbital elements')
disp(init_orbital)