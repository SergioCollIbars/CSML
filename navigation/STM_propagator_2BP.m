clear;
clc;
close all;
format long g;

addpath('functions/')
set(0,'defaultAxesFontSize',16);
%%          2BP STM PROPAGATOR
% Description: reconstruct trajectory based on ideal measurements and
% Newton Rapson integatrion method.

% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
tmin = 0;           % [s]
tmax = 1*86400;     % [s]
frec = 10/10;       % [Hz]
N = tmax*frec;
t = linspace(tmin, tmax, N);

% define planet parameters
R =  6.3781E6;  % [m]
m = 5.9722E24;  % [Kg]
GM = G*m;       % [m^3 s^-2]
n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
nmax = 0;
planetParams = [GM, R, nmax, 0];
poleParams = [2*pi/86400, deg2rad(180+45), deg2rad(164.11), deg2rad(89.374)];
Cmat = [1, 0, 0, 0;...
     0, 0, 0, 0;...
    -1.08E-3, 0, 1.57E-6, 0;...
     2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 

Smat = [0, 0, 0, 0;...
     0, 0, 0, 0;...
     0, 0, -9.03E-7, 0;...
     0, 2.68E-7, -2.12E-7, 1.98E-7]; 

% Initital conditions. Orbital parameters
e = 0.4;                     % eccentrycity
a = R + 1000e3;              % semi major axis [m]
rho = a * (1 - e^2);         % orbital param [m]

i = deg2rad(45);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(-180);           % true anomaly [rad]

[r01, v01] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
[r02, v02] = orbElems_2_ACI(rho, f+deg2rad(1), GM, Omega, omega, i, e);
[r03, v03] = orbElems_2_ACI(rho, f-deg2rad(10), GM, Omega, omega, i+deg2rad(0), e);

% define initial conditions. Small body
X01 = [r01(1); r01(2); r01(3); v01(1); v01(2); v01(3)];
X02 = [r02(1); r02(2); r02(3); v02(1); v02(2); v02(3)];
X03 = [r03(1); r03(2); r03(3); v03(1); v03(2); v03(3)];
X0 = X01;

% integatre trajectory and STM
STM_0 = reshape(eye(6,6), [36, 1]);

% define integration options
options = odeset('RelTol',1e-13,'AbsTol',1e-13);

% ODE 113
[time, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, "2BP"), t, [X0; STM_0], options);

% relative S/C propagation
[~, state_true1] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, "2BP"), t, [X01; STM_0], options);
[~, state_true2] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, "2BP"), t, [X02; STM_0], options);
[~, state_true3] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, "2BP"), t, [X03; STM_0], options);
STM1 = state_true1(:, 7:end);
STM2 = state_true2(:, 7:end);
STM = STM1;
STM3 = state_true3(:, 7:end);

% plot relative states
sigmar = 1E-10;  % [m]
sigmav = 1E-10;  % [m/s]
noiser = normrnd(0,sigmar, [3, N]);
noisev = normrnd(0,sigmav, [3, N]);
rho    = state_true2(:, 1:3)' - state_true1(:, 1:3)' + noiser; 
rhoDot = state_true2(:, 4:6)' - state_true1(:, 4:6)' + noisev;

rho1_RTN = ones(3, N) * NaN;
rho2_RTN = ones(3, N) * NaN;
beta_sc1 = ones(4, N) * NaN;
beta_sc2 = ones(4, N) * NaN;
lla_sc1 = ones(2, N) * NaN;
lla_sc2 = ones(2, N) * NaN;
for j  = 1:N
    NB_sc1 =  RTN2ECI(state_true1(j, 1:3)', state_true1(j, 4:6)');
    NB_sc2 =  RTN2ECI(state_true2(j, 1:3)', state_true2(j, 4:6)');

    T1 = NB_sc1(:, 2);
    T2 = NB_sc2(:, 2);

    BN_sc1 = NB_sc1';
    BN_sc2 = NB_sc2';
   
    rho1_RTN(:, j) = BN_sc1 * rho(:, j);
    rho2_RTN(:, j) = BN_sc2 * rho(:, j);

    beta_sc1(:, j) = DCM2EP(BN_sc1);
    beta_sc2(:, j) = DCM2EP(BN_sc2);
    
    p = sqrt(T1(1)^2 + T1(2)^2);
    phi = atan2(T1(3), p);
    lat = atan2(T1(2), T1(1));
    lla_sc1(:, j) = [phi; lat];

    p = sqrt(T2(1)^2 + T2(2)^2);
    phi = atan2(T2(3), p);
    lat = atan2(T2(2), T2(1));
    lla_sc2(:, j) = [phi; lat];
end

figure()
subplot(2, 1, 1)
plot(time./86400, vecnorm(rho)./1000, 'LineWidth', 2, 'Color', 'g')
ylabel('\rho [km]')
xlabel('TIME [days]')

subplot(2, 1, 2)
plot(time./86400, vecnorm(rhoDot)./1000, 'LineWidth', 2, 'Color', 'b')
ylabel('\rho dot [km/s]')
xlabel('TIME [days]')
sgtitle('relative measurements between S/C 1 and S/C 2')

% absolute state SC 1 and SC 2
figure()
subplot(2, 2, 1)
plot(time./86400, state_true1(:, 1:3)./1000, 'LineWidth', 2)
ylabel('r [km]')
xlabel('TIME [days]')
legend('r1', 'r2', 'r3')
title('SC 1')

subplot(2, 2, 3)
plot(time./86400, state_true1(:, 4:6)./1000, 'LineWidth', 2)
ylabel('v [km/s]')
xlabel('TIME [days]')
legend('v1', 'v2', 'v3')

subplot(2, 2, 2)
plot(time./86400, state_true2(:, 1:3)./1000, 'LineWidth', 2)
ylabel('r [km]')
xlabel('TIME [days]')
legend('r1', 'r2', 'r3')
title('SC 2')

subplot(2, 2, 4)
plot(time./86400, state_true2(:, 4:6)./1000, 'LineWidth', 2)
ylabel('v [km/s]')
xlabel('TIME [days]')
legend('v1', 'v2', 'v3')
sgtitle('relative measurements between S/C 1 and S/C 2')

% realtive motion RTN
figure()
for j = 1:3
    subplot(3, 1, j)
    plot(time./86400, rho1_RTN(j, :), 'LineWidth', 2)
    hold on;
    plot(time./86400, rho2_RTN(j, :), 'LineWidth', 2)
    xlabel('TIME [days]')
    ylabel('\rho ' + string(j))
    if(j == 1)
        legend('RNT SC1', 'RTN SC2')
    end
end
sgtitle('relative measurement in the RNT frame')

% reconstruct trajectory
gamma = zeros(9, N);
for j =1:N
    [~, ddU_ACI] = gradiometer_meas(t(j), state_true(j,1:3)', ...
    Cmat, Smat, poleParams, planetParams);
    gamma(:, j) = reshape(ddU_ACI, [9, 1]);
end

rint = zeros(3, N);
rint(:, 1) = state_true1(1, 1:3)';

rint_2 = zeros(3, N);
rint_2(:, 1) = state_true1(1, 1:3)';
vint = zeros(3, N);
w = vint;
for j = 1:N
    %vint(:, j) = trapz(state_true1(1:j, 4:6)', 2);
    %vint_2 = simps(state_true1(1:j, 4:6)', 2);
    %rint(:, j) = rint(:, 1) + vint(:, j);    
    %rint_2(:, j) = rint_2(:, 1) + vint_2;  
    w(:, j) = cross(state_true1(j, 1:3)', state_true1(j, 4:6)')./(vecnorm(state_true1(j, 1:3)')^2);
end

% propagate state
[state_nom] = propagateState(state_true, gamma, t);

% compute state derivative
Xdot_true = zeros(N, 6);
Xdot1_true = zeros(N, 6);
Xdot2_true = zeros(N, 6);

Xdot_true(:, 1:3) = state_true(:, 4:6);
Xdot1_true(:, 1:3) = state_true1(:, 4:6);
Xdot2_true(:, 1:3) = state_true2(:, 4:6);

for j = 1:N
    r = state_true(j, 1:3)';
    avec = - GM / (vecnorm(r)^3) * r;
    Xdot_true(j, 4:6) = avec';

    r1 = state_true1(j, 1:3)';
    avec = - GM / (vecnorm(r1)^3) * r1;
    Xdot1_true(j, 4:6) = avec';

    r2 = state_true2(j, 1:3)';
    avec = - GM / (vecnorm(r2)^3) * r2;
    Xdot2_true(j, 4:6) = avec';
end

% propagate the state derivative based on Hamiltonian dynamics
[state_nom, Xdot_nom] = propagateDynamics(Xdot1_true, state_true, gamma, t);

% relative states observability
% % O = zeros(3*N, 3);
% % H = [eye(3,3), zeros(3,3);zeros(3,3), eye(3,3)];
% % for  j = 1:N
% %     PHI1 = reshape(STM1(j, :), 6,6);
% %     PHI2 = reshape(STM2(j, :), 6,6);
% %     PHI3 = reshape(STM3(j, :), 6,6);
% %     %PHI_t = reshape(STM(j+1, :), 6,6);
% % 
% %     acc = Xdot_true(j, 4:6)';
% %     omega = w(:, j);
% %     [acc_tilde] = skewMatrix(omega);
% % 
% %     maxInd = 3*j;
% %     minInd  = maxInd - 2;
% % % %     O(minInd:maxInd, :) = [PHI1, -2*PHI1 + PHI2];
% % % %     O(minInd:maxInd, :) = [H * (PHI2 - PHI1), 0.5*H*(PHI1 + PHI2)];
% %    O(minInd:maxInd, :) = acc_tilde';
% % end
% % On = length(O(1, :));
% % l  = O'*O;
% % A  = l(1:On/2, 1:On/2);
% % B  = l(1:On/2, On/2+1:On);
% % C  = l(On/2+1:On, 1:On/2);
% % D  = l(On/2+1:On, On/2+1:On);
% % 
% % P11 = inv(A - B*inv(D)*C);
% % P12 = -P11*B*inv(D);
% % P21 = -inv(D)*C*P11;
% % P22 = inv(D) + inv(D)*C*P11*B*inv(D);
% % 
% % P = [P11, P12;P21, P22];
% % sigma = diag(P);

Xhat = zeros(6, N);
for j = 1:N
    PHI1 = reshape(STM1(j, :), 6,6);
    Xhat(:, j) = PHI1 *  [rho(:,1);rhoDot(:, 1)];
end

% compute initial dynamics
Ax = 0;
Nx = 0;
Rinv = [eye(3,3) * 1/(sigmar^2),zeros(3,3);...
    zeros(3,3), eye(3,3) * 1/(sigmav^2)];
At = time(2) - time(1);
for j =2:N-1
    % measurement
    y = [rho(:,j+1);rhoDot(:, j+1)] - [rho(:,j-1);rhoDot(:, j-1)];
    y = y./(2*At);
    %y = state_true2(j, 1:6)' - state_true1(j, 1:6)';
    %acc = Xdot1_true(j, 4:6)';
    %y = -cross(vint(:, j), acc);

    % visibility matrix
    PHI1 = reshape(STM1(j, :), 6,6);
    PHI2 = reshape(STM2(j, :), 6,6);
    H = [PHI2, -PHI1];
    %H = [PHI1,-2.*PHI1+PHI2];
    %H = [eye(6,6)*(PHI1 - PHI2), 0.5*eye(6,6)*(PHI1 + PHI2)];
    %H = -skewMatrix(acc);

    % y = PHI1 * X0 - y;
   
    Ax = Ax + H'* Rinv * H;
    Nx = Nx + H'* Rinv * y;
end
Xnot = Ax\Nx;
err = [Xdot2_true(1, :)';Xdot1_true(1, :)'] - Xnot;
P = diag(sqrt(inv(Ax)));

%% PLOT
% plot trajectory
t = time;
figure()
plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), ...
    'LineWidth', 2, 'Color', 'b')
axis equal
title('Real trajectory inertial frame')
grid on;
xlabel('[-]')
ylabel('[-]')
zlabel('[-]')

figure()
plot3(state_nom(:, 1), state_nom(:, 2), state_nom(:, 3), ...
    'LineWidth', 2, 'Color', 'g')
axis equal
title('Nominal trajectory inertial frame')
grid on;
xlabel('[-]')
ylabel('[-]')
zlabel('[-]')

% plot true and nominal states
ttitle1 = "true and nominal trajectory. Inertial frame";
ttitle2 = "absolute error state. Inertial frame";
plot_states(state_true, state_nom, t, ttitle1, ttitle2)

% plot true and nominal dynamics. 
ttitle1 = "true and nominal dynamics. Inertial frame";
ttitle2 = "absolute error dynamics. Inertial frame";
plot_states(Xdot_true, Xdot_nom, t, ttitle1, ttitle2)


%% functions
function [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e)
    % compute initial state vectors (orbit frame)
    r0 = rho / (1 + e * cos(f)) * [cos(f);...
                                  sin(f);...
                                  0];
    v0 = sqrt(GM / rho) * [-sin(f);...
                          e + cos(f);...
                          0];
    
    % compute rotation matrix: Body to ACI
    [BN] = rotationMatrix(Omega, i, omega, [3,1,3]);
    
    % rotate initial state vectors. {I, J, K} frame
    r0 = BN' * r0;
    v0 = BN' * v0;
end

function plot_states(state_true, state_nom, t, ttitle1, ttitle2)
    left = [1, 3, 5];
    right = [2, 4, 6];
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t, state_true(:, j),...
            'LineWidth', 2, 'Color', "#FF00FF")
        hold on;
        plot(t, state_nom(:, j),...
        'LineWidth', 2, 'Color', 'g')
    
        xlabel('TIME [-]')
        ylabel('r_' + string(j))
        if(j == 1)
            legend('true', 'recons')
        end
    
        subplot(3, 2, right(j));
        plot(t, state_true(:, j+3), ...
            'LineWidth', 2, 'Color', "#FF00FF")
        hold on;
        plot(t, state_nom(:, j+3), ...
            'LineWidth', 2, 'Color', 'g')
       xlabel('TIME [-]')
       ylabel('v_' + string(j))
    end
    sgtitle(ttitle1)
    
    figure()
    for j =1:3
        subplot(3, 2, left(j));
        plot(t./86400, abs(state_true(:, j) - state_nom(:, j)),...
            'LineWidth', 2, 'Color', "r")
        xlabel('TIME [-]')
        ylabel('r_' + string(j))
        legend('true - recons')
    
        subplot(3, 2, right(j));
        plot(t./86400, abs(state_true(:, j+3) - state_nom(:, j+3)), ...
            'LineWidth', 2, 'Color', "r")
       xlabel('TIME [-]')
       ylabel('v_' + string(j))
       legend('true - recons')
    end
    sgtitle(ttitle2)
end

%%              FUNCTIONS

function [state_nom] = propagateState(state_true, gamma, t)
    deltaX0 = (state_true(2, 1:6) - state_true(1, 1:6))';
    [X_new] = reconstruct_traj(t, gamma, deltaX0, "2BP", 4);
    
    state_nom = state_true * 0;
    state_nom(1, 1:6) = state_true(1, 1:6);
    for j =2:length(t)
        %PHI = reshape(state_true(j, 7:end), [6,6]);
        %deltaX = PHI * deltaX0;
        deltaX = X_new(:, j);
        state_nom(j, 1:6) = state_nom(j-1, 1:6) + deltaX';
    end
end

function [state_nom, Xdot_nom] = propagateDynamics(Xdot_true, state_true, gamma, t)
    N = length(t);
    % compute Dynamics from STM and initial value [v0; a0]
    Xdot_nom = Xdot_true * 0;
    Xdot_nom(1, :) = Xdot_true(1, :);
    for j = 2:N
        PHI = reshape(state_true(j, 7:end), [6,6]);
        R1 = PHI;
        Xdot_nom(j, :) = (R1 * Xdot_nom(1, :)')';
    end

    % reconstruc state from dynamics. Taylor expansion
    state_nom = zeros(N, 6);
    state_nom(1, :) = state_true(1, 1:6);
    At = t(2) - t(1);
    for j = 1:N-1
        T = reshape(gamma(:, j), [3,3]);
        J = [zeros(3,3), eye(3,3);T,zeros(3,3)];
        dX = Xdot_nom(j, :)';
        ddX = J * dX;
        if(j>1)
            Tdot = (reshape(gamma(:, j+1), [3,3]) - reshape(gamma(:, j-1), [3,3]))./(2*At); 
            Jdot = [zeros(3,3), zeros(3,3);Tdot,zeros(3,3)];
            dddX = Jdot*dX + J*ddX;
        else
            dddX = 0;
        end
        
        state_nom(j+1, :) = state_nom(j, :) + (dX')*At + 1/2*(ddX')*At*At + ...
            1/6*(dddX')*At*At*At;
    end
end

function [r] = radialVec_solver(GM, gamma)
    N = length(gamma(1, :));
    r = zeros(3, N);
    for j = 1:N
        T = reshape(gamma(:, j), [3,3]);
        [V, D] = eig(T);
        [e, p] = max(max(D));    % max eigenvalue
        n = V(:, p);

        r(:, j) = (2*GM/ e)^(1/3) * abs(n);
        rx = r(1, j);
        ry = r(2 ,j);
    end
end

function [x_tilde] = skewMatrix(x)
    x_tilde = [0, -x(3), x(2);...
            x(3), 0, -x(1);...
            -x(2), x(1), 0];
end

function [RTN_2_ECI] = RTN2ECI(r, v)
    % ------------------------------------------------------------------- %
    %                     RTN 2 ECI FRAME FUNCTION
    % Author: Sergio Coll Ibars

    % Date: 02/09/2023

    % Description: Compute the radial, along-track, cross-track (RTN) 
    %   to Earth-centered inertial rotation matrix. If applied to a
    %   position vector in the RTN frame, it will transform that vector to
    %   into the equivalent position vector in the ECI frame.

    % Input:
    %   r: current position vector. ECI frame
    %   v: current velocity vector. ECI frame

    % Output:
    %   RTN_2_ECI: rotation matrix from RTN to ECI frame
    % ------------------------------------------------------------------- %

    n = cross(r, v);

    R = r / vecnorm(r);
    N = n / vecnorm(n);
    T = cross(N, R);

    RTN_2_ECI = [R, T, N];

end

function [beta] = DCM2EP(C)
    b0sq = 0.25*(1 + trace(C));
    b1sq = 0.25*(1 + 2*C(1,1) - trace(C));
    b2sq = 0.25*(1 + 2*C(2,2) - trace(C));
    b3sq = 0.25*(1 + 2*C(3,3) - trace(C));

    X = [b0sq; b1sq; b2sq; b3sq];
    [~, ind] = max(X);

    if(ind == 1)        % b0 max
        b0 = sqrt(b0sq);
        b1 = (C(2,3) - C(3,2))/(4*b0);
        b2 = (C(3,1) - C(1,3))/(4*b0);
        b3 = (C(1,2) - C(2,1))/(4*b0);
    elseif(ind == 2)    % b1 max
        b1 = sqrt(b1sq);
        b0 = (C(2,3) - C(3,2))/(4*b1);
        b2 = (C(1,2) - C(2,1))/(4*b1);
        b3 = (C(3,1) - C(1,3))/(4*b1);
    elseif(ind == 3)    % b2 max
        b2 = sqrt(b2sq);
        b0 = (C(3,1) - C(1,3))/(4*b2);
        b1 = (C(1,2) - C(2,1))/(4*b2);
        b3 = (C(2,3) - C(3,2))/(4*b2);
    elseif(ind == 4)    % b3 max
        b3 = sqrt(b3sq);
        b0 = (C(1,2) - C(2,1))/(4*b3);
        b1 = (C(3,1) - C(1,3))/(4*b3);
        b2 = (C(2,3) - C(3,2))/(4*b3);
    end
    beta = [b0;b1;b2;b3];
end