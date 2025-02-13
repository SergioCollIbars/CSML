clc;
clear;
close all;
format long g;

addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);
%%               ATTITUDE ERROR ANALYSIS V2
%   Description: Analysis on how the attitude influence the measurements.
%   test under RTN frame or eigenvectors frame.

Al = sym("Al","real");
Ap = sym("Ap","real");
Ar = sym("Ar","real");
G1 = sym("G1","real");
G2 = sym("G2","real");
G3 = sym("G3","real");

G11 = sym("G11","real");
G12 = sym("G12","real");
G13 = sym("G13","real");

G21 = sym("G21","real");
G23 = sym("G23","real");

M3 = [cos(Al), sin(Al), 0;-sin(Al), cos(Al), 0;0, 0, 1];
M2 = [cos(Ap), 0, -sin(Ap); 0, 1, 0;sin(Ap), 0, cos(Ap)];

M1 = [1, 0, 0;0, cos(Ar), sin(Ar);0, -sin(Ar), cos(Ar)];
%M1 = [1,0,0;0,cos(pi),0;0,0,cos(pi)];
C = M1 * M2 ;
%C = M1 * M2 * M3;
C = M3;

T =[G1,0,0;0,G2,0;0,0,G3];
T = C' * T * C;

% Harmonic values. A priori
Cnm = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0.1 0.2 0 0 0 0;...
        0.4 0 0 0 0 0 0;...
        0.02 0.0001 -0.0004 -7-05 4.51-05 0 0;...
        -0.005 0.00044 9-05 1e-05 1e-05 1e-06 0;...
        -0.009 0.001 9e-05 2e-05 -2e-06 9e-08 2e-07];
 
Cnm = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0 0 0 0 0 0;...
        0.4 0 0 0 0 0 0;...
        0.02 0 0 0 0 0 0;...
        -0.005 0 0 0 0 0 0;...
        -0.009 0 0 0 0 0 0];

Snm = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;
        0 0.007 0.0001 0.0001 0.003 0 0;...
        0 -2e-05 -8e-05 5e-06 -0.0005 0.0005 0;...
        0 -1e-06 -1e-05 -8e-06 -2e-06 -0.0006 -4e-05];
Snm = zeros(7,7);

% Bennu parameters
GM = 5.2;
Re = 246;
n = sqrt(GM / 1000^3);
n_max = 6;
W = 4.06130329511851E-4;
W0 = 0;
RA = deg2rad(86.6388);
DEC = deg2rad(-65.1086);

% number of points
t0 = 0;
t_max = round(0.5*86400);
Nt_max = round(t_max / 1);
t = linspace(t0, t_max, Nt_max);

% Inital conditions
e = 0;                       % eccentrycity
a = 500;             % semi major axis [m]
rho = a * (1 - e^2);         % orbital param [m]

i = deg2rad(45);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(0);           % true anomaly [rad]

[r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
X0 = [r0;v0];
dX = [1;1;1;0;0;0];

% integrate nominal and true state
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, 0, W0, W, ...
    RA, DEC), t, X0, options);
[~, state_n] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, 0, W0, W, ...
    RA, DEC), t, X0 + dX, options);


% compute measurements in both orbits
Y_t = zeros(9, Nt_max);
Y_n = Y_t;
ang1 = zeros(3, Nt_max);    % truth EIG frame angle
ang2 = zeros(3, Nt_max);    % nominal EIG frame angle
ang3 = zeros(3, Nt_max);    % truth RTN frame angle
for k = 1:Nt_max
    rt_ACI = state_t(k, 1:3)';
    vt_ACI = state_t(k, 4:6)';
    rn_ACI = state_t(k, 1:3)' + [10;10;10];
    %rn_ACI = state_n(k, 1:3)';
    vn_ACI = state_n(k, 4:6)';

    [ACI_RTN] = RTN2ECI(rt_ACI, vt_ACI);
    
    Wt = W0 + W * t(k);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
    r_ACAF = ACAF_ACI * rt_ACI;
    ru_ACAF = r_ACAF./vecnorm(r_ACAF);
    [ACAF_ENU] = ENU2ACAF(ru_ACAF);

    % ACAF coordinates
    rt_ACAF = ACAF_ACI * rt_ACI;
    rn_ACAF = ACAF_ACI * rt_ACI + [10;10;10];
    
    % compute truth meas. rotated frame frame
    [~, ~, ddU_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rt_ACAF, Re, GM, 0);
    ddU_ACI = ACAF_ACI' * ddU_ACAF * ACAF_ACI; % ACI
    ddU = ACI_RTN'  * ddU_ACI      * ACI_RTN; % RTN
    %ddU = ACAF_ENU' * ddU_ACAF * ACAF_ENU;    % ENU
    ddU = ddU_ACAF;
    %ddU = ddU_ACI;

    [ev, v] = eig(ddU);
    L = ev;
    %L = eye(3,3);
    
% %     rv_ACAF = L(:, 3);
% %     v_ACAF = ACAF_ACI * vn_ACI;
% %     [ACAF_RTN] = RTN2ECI(rv_ACAF, v_ACAF);
% %     L = ACAF_RTN;

    % save angles
    ev = ev';
    ang1(1, k) = atan2(ev(1, 2), ev(1,1));
    ang1(2, k) = -asin(ev(1, 3));
    ang1(3, k) = atan2(ev(2, 3), ev(3,3));

    ang3(1, k) = atan2(ACAF_ENU(1, 2), ACAF_ENU(1,1));
    ang3(2, k) = - asin(ACAF_ENU(1, 3));
    ang3(3, k) = atan2(ACAF_ENU(2, 3), ACAF_ENU(3,3));

% %     ang1(:, k) = ang3(:, k);

    Y_t(:, k) = reshape(L' * ddU * L, [9, 1]);

    % compute nominal meas. rotated reference frame
    [~, ~, ddU_ACAF] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rn_ACAF, Re, GM, 0);
    ddU_ACI = ACAF_ACI' * ddU_ACAF * ACAF_ACI; % ACI
    ddU = ACI_RTN'  * ddU_ACI * ACI_RTN;      % RNT
    %ddU = ACAF_ENU' * ddU_ACAF * ACAF_ENU;    % ENU
    ddU = ddU_ACAF;
    %ddU = ddU_ACI;
    [ev2, v2] = eig(ddU);

    Y_n(:, k) = reshape(L' * ddU * L, [9, 1]);
    ev2 = ev2';

    ang2(1, k) = atan2(ev2(1, 2), ev2(1,1));
    ang2(2, k) = - asin(ev2(1, 3));
    ang2(3, k) = atan2(ev2(2, 3), ev2(3,3));
     b = abs(ang2(3, k) - ang1(3, k));
    if ((deg2rad(.1) < b)  && (b < deg2rad(1)))
        [partials] = computePartials(rn_ACI, Cnm, Snm, GM, Re, n_max, ...
            ev, ACAF_ACI);
        %%disp(partials)
    else
% %         [partials] = computePartials(rn_ACI, Cnm, Snm, GM, Re, n_max, ...
% %             ev, ACAF_ACI);
% %         disp(partials)
    end
end
err = Y_t - Y_n;

ang1 = wrapTo2Pi(ang1);
ang2 = wrapTo2Pi(ang2);
ang3 = wrapTo2Pi(ang3);

figure()
lw = 2;
for k = 1:3
    subplot(3, 1 ,k);
    plot(t, rad2deg(ang1(k, :)), t, rad2deg(ang2(k, :)), 'LineWidth', lw)
end
sgtitle('EIG frame attitude comparison')
legend('truth', 'nominal')

figure()
for k = 1:3
    subplot(3, 1 ,k);
    plot(t, rad2deg(ang1(k, :) - ang2(k, :)), 'LineWidth', lw)
end
sgtitle('Error EIG frame attitude')

figure()
for k = 1:3
    subplot(3, 1 ,k);
    plot(t, rad2deg(ang3(k, :)), 'LineWidth', lw)
end
sgtitle('Truth RTN frame attitude')

% plot orbit
plotOrbit(state_t(:, 1:3)', state_n(:, 1:3)', Nt_max);

% plot measurements
ttitle = 'error along orbit. EIG frame';
ytitle = ["R", "RT", "RN", "T", "TN", "N"];
plotT(err, Nt_max, t, ttitle, ytitle)

ttitle = 'truth meas. EIG frame';
ytitle = ["R", "RT", "RN", "T", "TN", "N"];
plotT(Y_t, Nt_max, t, ttitle, ytitle)

ttitle = 'nominal meas. truth EIG frame';
ytitle = ["R", "RT", "RN", "T", "TN", "N"];
plotT(Y_n, Nt_max, t, ttitle, ytitle)

% plot position error
figure()
for k = 1:3
    subplot(3, 1, k)
    plot(t, state_t(:, k)' - state_n(:, k)', 'LineWidth', 2, ...
        'Color', 'g')
    ylabel('error ' + string(k))
end
sgtitle('Position error in meters')

%%  FUNCTIONS
function plotOrbit(r_ACI, rP_ACI, Nt)
    % plot 3D orbit
    figure();
    f = zeros(1, 2);
    f(1) = plot3(r_ACI(1,1:Nt), r_ACI(2,1:Nt), r_ACI(3,1:Nt), 'MarkerSize', 10, ...
        'LineWidth', 2);
    hold on;
    f(2) = plot3(rP_ACI(1, 1:Nt), rP_ACI(2, 1:Nt), rP_ACI(3, 1:Nt), ...
        'MarkerSize', 10, 'LineWidth', 2);
    legend(f, 'Nominal', 'Perturbed')
    title('Orbit in the inertial frame {I, J, K}');
    grid on;
    axis equal;
    
    % plot axis
    mAxis = max(max(r_ACI));
    axis([0 mAxis 0 mAxis 0 mAxis])
    hold all;
    quiver3(0,0,-max(0),0,0,max(zlim),'r','LineWidth',1)
    quiver3(0,-max(0),0,0,max(ylim),0,'r','LineWidth',1)
    quiver3(-max(0),0,0,max(xlim),0,0,'r','LineWidth',1)
    text(0,0,max(zlim),'K','Color','r')
    text(0,max(ylim),0,'J','Color','r')
    text(max(xlim),0,0,'I','Color','r')
    
    scale =  450;  % Bennu object scale factor
    obj = readObj('Bennu-Radar.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 
    
    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
end


function plotT(Err, Nt_max, TIME, ttitle, ytitle)
    % time vector [h]
    t = TIME(1:Nt_max)./ 87132.103;
    lw = 2;
    cl = 'r';

    figure();

    subplot(2, 3, 1);
    plot(t, Err(1, :), 'LineWidth', lw, 'Color', cl)
    xlabel("Orbital rev, T")
    ylabel(ytitle(1))
    grid on;

    subplot(2, 3, 2);
    plot(t, Err(2, :), 'LineWidth', lw, 'Color', cl)
    xlabel("Orbital rev, T")
    ylabel(ytitle(2))
    grid on;

    subplot(2, 3, 3);
    plot(t, Err(3, :), 'LineWidth', lw, 'Color', cl)
    xlabel("Orbital rev, T")
    ylabel(ytitle(3))
    grid on;

    subplot(2, 3, 4);
    plot(t, Err(5, :), 'LineWidth', lw, 'Color', cl)
    xlabel("Orbital rev, T")
    ylabel(ytitle(4))
    grid on;

    subplot(2, 3, 5);
    plot(t, Err(6, :), 'LineWidth', lw, 'Color', cl)
    xlabel("Orbital rev, T")
    ylabel(ytitle(5))
    grid on;

    subplot(2, 3, 6);
    plot(t, Err(9, :), 'LineWidth', lw, 'Color', cl)
    xlabel("Orbital rev, T")
    ylabel(ytitle(6))
    grid on;

    sgtitle(ttitle)
end

function [ACAF_ENU] = ENU2ACAF(r_ACAF)
    rxy = sqrt(r_ACAF(1)^2 + r_ACAF(2)^2);

    phi = atan2(r_ACAF(3), rxy);
    l = atan2(r_ACAF(2), r_ACAF(1));

    ACAF_ENU = [-sin(l), -sin(phi)*cos(l), cos(phi)*cos(l);...
        cos(l), -sin(phi)*sin(l), cos(phi)*sin(l);...
        0, cos(phi), sin(phi)];
end

function [partials] = computePartials(r0, Cnm, Snm, GM, Re, n_max, ...
    EIG_ACAF, ACAF_ACI)
    % output matrix
    partials = zeros(9, 3);

    % set correctly the coefficient matrix
    Cnm(:, 2:end) = Cnm(:, 2:end).*0; 

    % compute initial spherical coordinates
    r0 = ACAF_ACI * r0; % ACAF coordiantes
    r0n = r0./vecnorm(r0);
    
    phi0 = pi/2 - atan2(r0n(3), sqrt(r0n(1)^2 + r0n(2)^2));
    lambda0 = atan2(r0n(2), r0n(1));
    rho0 = vecnorm(r0);
    
    Asph = zeros(3, 1);
    for j = 1:3
        % compute displaced spherical coordinates
        Asph(j) = 1E-6;
        rho    = rho0 + Asph(1);
        phi    = phi0 + Asph(2);
        lambda = lambda0 + Asph(3);
        
        rx = rho * cos(phi) * cos(lambda);
        ry = rho * cos(phi) * sin(lambda);
        rz = rho * sin(phi);
        
        rpos = [rx;ry;rz];  % positive increment, ACAF
        
        rho    = rho0 - Asph(1);
        phi    = phi0 - Asph(2);
        lambda = lambda0 - Asph(3);
        
        rx = rho * sin(phi) * cos(lambda);
        ry = rho * sin(phi) * sin(lambda);
        rz = rho * cos(phi);
        
        rneg = [rx;ry;rz];  % negative increment, ACAF

        [~, ~, ddU_ACAF_pos] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                           rpos, Re, GM, 0);


        [~, ~, ddU_ACAF_neg] = potentialGradient_nm(Cnm, Snm, n_max, ...
                                                   rneg, Re, GM, 0);

        % project into the EIG frame
        ddU_EIG_pos = EIG_ACAF * ddU_ACAF_pos * EIG_ACAF';
        ddU_EIG_neg = EIG_ACAF * ddU_ACAF_neg * EIG_ACAF';

        rpos = EIG_ACAF * rpos;
        rneg = EIG_ACAF * rneg;

        y_pos = reshape(ddU_EIG_pos, [9, 1]);
        y_neg = reshape(ddU_EIG_neg, [9, 1]);
        partials(:, j) = (y_pos - y_neg)./(vecnorm(rpos-rneg));
    end
end