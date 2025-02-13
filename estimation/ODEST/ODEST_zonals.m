clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);
%%              ODEST ALGORTIHM INCLUDING ZONALS
% Description: modify the ODEST algorithm to work eliminate error in the
% zonals coefficients too.


% gravity field coefficients
Cnm_truth = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0.1 0.2 0 0 0 0;...
        0.2 0.3 0.04 0 0 0 0;...
        0.02 0.0001 -0.0004 -7-05 4.51-05 0 0;...
        -0.005 0.00044 9-05 1e-05 1e-05 1e-06 0;...
        -0.009 0.001 9e-05 2e-05 -2e-06 9e-08 2e-07];

Snm_truth = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 1E-6 -3e-05 0 0 0 0;...
        0 0.0001 7e-05 -0.001 0 0 0;
        0 0.007 0.0001 0.0001 0.003 0 0;...
        0 -2e-05 -8e-05 5e-06 -0.0005 0.0005 0;...
        0 -1e-06 -1e-05 -8e-06 -2e-06 -0.0006 -4e-05];

Cnm_nom = Cnm_truth;
Cnm_nom(:, 2:end) = Cnm_nom(:, 2:end) * 0;
Snm_nom = zeros(7,7);

Cnm_nom = Cnm_truth;
Snm_nom = Snm_truth;

Cnm_dist = Cnm_truth-Cnm_nom;
Snm_dist = Snm_truth-Snm_nom;

% Bennu parameters
GM = 5.2;
Re = 246;
n_max = 4;

% number of points
t0 = 0;
t_max = round(10*86400);

% Inital conditions
e = 0;                       % eccentrycity
a = 2000;             % semi major axis [m]
rho = a * (1 - e^2);         % orbital param [m]
W0 = 0;
W = 0;
RA = -pi/2;
DEC = pi/2;

i = deg2rad(40);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(0);           % true anomaly [rad]

[r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
X0 = [r0;v0];

% integrate nominal and true state
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[t, state_t] = ode113(@(t, x) EoM(t, x, Cnm_truth, Snm_truth, n_max, GM, Re, 0, W0, W, ...
    RA, DEC), [t0, t_max], X0, options);
Nt_max = length(t);

% compute measurement difference in both orbits.
Y_RTN = ones(6, Nt_max)* NaN;
distY = ones(6, Nt_max) * NaN;
Ar_ECEF = [1;1;1]*1E-2;
for k = 1:Nt_max
    % truth & nominal position 
    rt_ECEF = state_t(k, 1:3)';
    vt_ECEF = state_t(k, 4:6)';
    rn_ECEF = rt_ECEF + Ar_ECEF;

    % truth and nominal meas
    [~, ~, Tt_ECEF] = potentialGradient_nm(Cnm_truth, Snm_truth, n_max, ...
                                           rt_ECEF, Re, GM, 0);
    [~, ~, Tdt_ECEF] = potentialGradient_nm(Cnm_dist, Snm_dist, n_max, ...
                                           rt_ECEF, Re, GM, 0);
    [~, ~, Tdn_ECEF] = potentialGradient_nm(Cnm_nom, Snm_nom, n_max, ...
                                           rn_ECEF, Re, GM, 0);
    
    % compute RTN frame
    [ECEF_RTN] = RTN2ECI(rn_ECEF, vt_ECEF);
    rn_RTN = ECEF_RTN' * rn_ECEF;
    Ar_RTN = ECEF_RTN' * Ar_ECEF; 

    % compute truth and nominal meas in nominal RTN frame
    Tt_RTN = ECEF_RTN' * Tt_ECEF * ECEF_RTN;
    Tn_RTN = ECEF_RTN' * Tdn_ECEF * ECEF_RTN;

    % save values
    dY_RTN = Tt_RTN - Tn_RTN;

    Y_RTN(1, k) = dY_RTN(1, 1);
    Y_RTN(2, k) = dY_RTN(1, 2);
    Y_RTN(3, k) = dY_RTN(1, 3);
    Y_RTN(4, k) = dY_RTN(2, 2);
    Y_RTN(5, k) = dY_RTN(2, 3);
    Y_RTN(6, k) = dY_RTN(3, 3);  

    % compute position partial. 1st order
    [Hpos] = compute_posPartials(n_max, 0, Cnm_nom, Snm_nom, Re, GM, rn_RTN);

    % compute Point mass approx error
    [Hpos0] = compute_posPartials(0, 0, Cnm_nom, Snm_nom, Re, GM, rn_RTN);
    Hposc = Hpos - Hpos0;
    distY(:, k) = (Hposc * Ar_RTN);
end

% plots
lab =  {'\Gamma_{RR}', '\Gamma_{RT}', '\Gamma_{RN}',...
        '\Gamma_{TT}', '\Gamma_{TN}', '\Gamma_{NN}'};
plot_QGG_components(t, Y_RTN, "Gradiometer components error. RTN frame", ...
   lab);

plot_trajectory(state_t);

figure()
plot(t./86400, distY, 'lineStyle', '-', 'LineWidth', 2)
title('Point mass approximation error')
legend(lab)

%% FUNCTIONS
function [ECEF_ENU] = rotate2ENU(p, l)
    ECEF_ENU = [-sin(l), -sin(p)*cos(l), cos(p)*cos(l);...
        cos(l), -sin(p)*sin(l), cos(p)*sin(l);...
        0, cos(p), sin(p)];
end

function [T]  = compute_GGT(r, p, l, GM, J2, J3, R) 
    % Description: return the grav. gradient tensor compoenents in
    % ENU coordinates including the J2 and J3 effect
    r3 = r^3;
    r5 = r^5;
    r6 = r^6;
    A = J2 * GM * R^2;
    B = J3 * GM * R^3;

    Trr = 2*GM/r3 - 6 * A / r5 * (3*sin(p)^2 - 1) + ...
        -10 * B / r6 * (5*sin(p)^3 - 3*sin(p));
    Trp = 6 * A /r5 * sin(2*p) + ...
        5/2 * B /r6 * (15*sin(p)*sin(2*p)/2 - 3*cos(p));
    Trl = 0;
    
    Tpp = -GM/r3 +3/2 * A /r5 * (7*sin(p)^2 - 3) + ...
        -1/2 * B /r6 * (3*sin(p) + 15*cos(p)*sin(2*p)/2 +30*sin(p)*cos(2*p)/2) + ...
        2 * B /r6 * (5*sin(p)^3 - 3*sin(p));
    Tpl = 0;

    Tll = -(Trr + Tpp);
    T = [Tll ,Tpl, Trl;Tpl, Tpp, Trp;Trl, Trp, Trr];
end

function [] = plot_QGG_components(t, Y, ttitle, llegend)
    figure()
    for j = 1:6
        subplot(2, 3, j)
        plot(t./86400, Y(j, :), 'LineWidth', 2)
        xlabel('TIME  [days]')
        ylabel('[1/s^2]')
        title(llegend(j))
    end
    sgtitle(ttitle)
end

function [] = plot_trajectory(state_t)
    figure()
    plot3(state_t(:,1), state_t(:,2), state_t(:,3), 'LineWidth', 2, ...
        'Color','r')
    axis equal;
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('Z [m]')
    title('S/C trajectory')
    grid on;
    % plot axis
    mAxis = max(max(state_t(:, 1:3)));
    axis([0 mAxis 0 mAxis 0 mAxis])
    hold all;
    quiver3(0,0,-max(0),0,0,max(zlim),'b','LineWidth',1)
    quiver3(0,-max(0),0,0,max(ylim),0,'b','LineWidth',1)
    quiver3(-max(0),0,0,max(xlim),0,0,'b','LineWidth',1)
    text(0,0,max(zlim),'K','Color','b')
    text(0,max(ylim),0,'J','Color','b')
    text(max(xlim),0,0,'I','Color','b')
    
    scale =  450;  % Bennu object scale factor
    obj = readObj('Bennu-Radar.obj');
    p = obj.v * 2 * scale;
    f = obj.f.v ; 
    
    trisurf(f,p(:,1),p(:,2),p(:,3));
    colormap(gray);
    axis equal
end