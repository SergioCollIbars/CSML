clear;
clc;
close all;

addpath('functions/')
set(0,'defaultAxesFontSize',16);
%%              STM PROPAGATION, COMPUTE COST FUNCTION
% Description: Look for minima in the cost function numerically

% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
tmin = 0;           % [s]
tmax = 20*500;   % [s]
frec = 10/10;        % [Hz]
N = tmax*frec;
t = linspace(tmin, tmax, N);

% define planet parameters
R =  6.3781E6;  % [m]
m = 5.9722E24;  % [Kg]
GM = G*m;       % [m^3 s^-2]
n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
nmax = 3;
planetParams = [GM, R, nmax, 0];
poleParams = [2*pi/86400, 0, -pi/2, pi/2];
Cmat = [1, 0, 0, 0;...
     0, 0, 0, 0;...
    -1.08E-3, 0, 1.57E-6, 0;...
     2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 

Smat = [0, 0, 0, 0;...
     0, 0, 0, 0;...
     0, 0, -9.03E-7, 0;...
     0, 2.68E-7, -2.12E-7, 1.98E-7]; 

% Initital conditions. Orbital parameters
e = 0.6;                     % eccentrycity
a = R + 3e3;             % semi major axis [m]
rho = a * (1 - e^2);       % orbital param [m]

i = deg2rad(90);            % inclination [rad]
omega = deg2rad(0);        % arg periapsis [rad]
Omega = deg2rad(0);        % RAAN [rad]
f = deg2rad(-180);         % true anomaly [rad]

[r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);

% define initial conditions. Small body
X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
deltaX = [100;0;0;0;0;0];
STM_0 = reshape(eye(6,6), [36, 1]);

% define integration options
options = odeset('RelTol',1e-13,'AbsTol',1e-13);

% ODE 113
[time, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, "2BP"), t, [X0; STM_0], options);

[~, state_true2] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
    poleParams, Cmat, Smat, "2BP"), t, [X0+ deltaX; STM_0], options);

% compute dynamics vector
[stateDot_true] = computeDynamic_vector(time, state_true, planetParams,...
    poleParams, Cmat, Smat);

[stateDot_true2] = computeDynamic_vector(time, state_true2, planetParams,...
    poleParams, Cmat, Smat);

% plot orbit states
plotOrbit(time, state_true, stateDot_true);

% compute cost function
N = 20;
J = zeros(N, N);
Jx = zeros(N, N, 6);
Jy = zeros(N, N, 6);

xAxis = zeros(1, N/4);
yAxis = zeros(1, N);
STM = state_true(:, 7:end);
deltaX = [1;1;1;1;1;1];
s = 2;
d = length(time);
X0t = stateDot_true(s, 1:6)';
X1t = stateDot_true(s+1, 1:6)';
[X, Y] = meshgrid(1:N, 1:N);
x0g = normrnd(0,1,[6,N]);
x1g = normrnd(0,1,[6,N]);
for i = 1:N % x direction
    for j = 1:N % y direction
        P = reshape(STM(s, :), [6,6]); 
        PP = reshape(STM(d, :), [6,6]); 
        p11 = P(1:3,1:3);
        p12 = P(1:3,4:6);
        p21 = P(4:6,1:3);
        p22 = P(4:6,4:6);

        B = P'*P; 
        b11 = B(1:3,1:3);
        b12 = B(1:3,4:6);
        b21 = B(4:6,1:3);
        b22 = B(4:6,4:6);

% %         X0 = stateDot_true(1, 1:6)' + 1*(i - N/2) * [0;0;0;deltaX(4:6)];
        X0 = stateDot_true(1, 1:6)' + 1*(i - N/2) * deltaX;
        X1 = stateDot_true(s, 1:6)' + 1*(j - N/2) * deltaX;

% %         X0 = stateDot_true(1, 1:6)' + 1*(i - N/2) * x0g(:, i);
% %         X1 = stateDot_true(s, 1:6)' + 1*(j - N/2) * x1g(:, j);

        %E = X0'*X1 - X0'*(P * X0);
        E = X1 - (P * X0);
        
        % plot axis
        xsign = (i-N/2) / abs(i-N/2);
        ysign = (j-N/2) / abs(j-N/2);
        xAxis(i) = xsign * round(vecnorm(X0 - stateDot_true(1, 1:6)'));
        yAxis(j) = ysign * round(vecnorm(X1 - stateDot_true(s, 1:6)'));

        % cost funciton value
        J(j, i)  = 0.5 * (E'*E);

        % function gradient
        A = P'*P;
        gx = (-P'*X1 + 0.5*(A + A')*X0);
        gy = (X1 - P*X0);
        for t = 1:6
            Jx(j, i, t) = gx(t);
            Jy(j, i, t) = gy(t);
        end
        
        % constuct Hesian matrix
       
        mu = vecnorm([gx;gy]);
        H = [0.5*(A+A'), -P';-P, eye(6,6)];
        G = H + mu * eye(12,12);
        f = H - H';
        D = eig(H);
        
        Hh = [0.5*(PP'*P + P'*PP), -0.5*PP', -0.5*P';...
            -0.5*PP, zeros(6,6), 0.5*eye(6,6);...
            -0.5*P, 0.5*eye(6,6), zeros(6,6)];
       
        % assume v0 known
        Hv0 = [b22, -p12, -p22;-p12, eye(3,3), zeros(3,3);...
            -p22, zeros(3,3), eye(3,3)];
        fv0 = Hv0 - Hv0';
        Dv0 = eig(Hv0);

        % assume v0, v1 known
        Hv = [b22,-p22;-p22,eye(3,3)];
        fv = Hv'- Hv;
        Dv = eig(Hv);
    end
end

% minumum function
minimum = min(min(J));
[Ix,Iy]=find(J==minimum);

% plot cost function
figure()
surf(X, Y, J)
xlabel("|\delta X_0|")
ylabel("|\delta X_i|")
%zlim([0, 5E3])
title(' Cost function mapping X0 -> Xi, i = ' + string(s));
%shading interp
set(gca, 'XTick',1:N, 'XTickLabel',...
    xAxis)
set(gca, 'YTick',1:N, 'YTickLabel',...
    yAxis)
hold on
plot3(X(N/2, N/2), Y(N/2, N/2), J(N/2, N/2), 'Marker','x', 'Color','k', ...
    'LineWidth', 4)
hold on;
plot3(X(Ix, Iy), Y(Ix, Iy), J(Ix, Iy), 'Marker','x', 'Color','r', ...
    'LineWidth', 4)
legend('cost function, J', 'true state', 'min(J) = ' + string(minimum))

% plot gradient
figure()
cx = zeros(N, N);
cy = cx;
for k =1:6
    cx = Jx(:, :, k);
    cy = Jy(:, :, k);
    subplot(2, 3, k)
    quiver(X, Y, cx, cy, 'Color', 'r', "LineWidth", 2)
end

%%      FUNCTIONS
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

function [Xdot] = computeDynamic_vector(time, state_true, planetParams,...
    poleParams, C_mat, S_mat)
    N  = length(time);
    Xdot = zeros(N, 6);

    % planet varaibles
    GM = planetParams(1);
    Re = planetParams(2);
    n_max = planetParams(3);
    Normalized = planetParams(4);

    W = poleParams(1);
    W0 = poleParams(2);
    RA = poleParams(3);
    DEC = poleParams(4);
    
    % velocity state
    Xdot(:, 1:3) = state_true(:, 4:6);
    
    % acceleration state
    for j = 1:N
        % ACI position vector
        r_ACI = [state_true(j, 1); state_true(j, 2); state_true(j, 3)];

        % Inertial to ACAF rotation 
        Wt = W * time(j) + W0;
        ACAF_N = rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % ACAF position vector
        r_ACAF = ACAF_N * r_ACI;

        [~, dU, ~] = potentialGradient_nm(C_mat, S_mat, n_max, ...
                                                    r_ACAF, Re, GM, ...
                                                    Normalized);

        % rotate from ECEF 2 ECI
        dU = ACAF_N' * dU;

        Xdot(j, 4:6) = dU'; 
    end
end

function plotOrbit(time, state_true, stateDot_true)
    figure()
    plot3(state_true(:, 1), state_true(:, 2), state_true(:, 3), ...
        'LineWidth', 2, 'Color', 'r');
    title('Orbit. Inertial frame')
    axis equal
    grid on;
    
    figure()
    for j = 1:3
        subplot(3, 1, j)
        plot(time, state_true(:, j), 'LineWidth', 2, 'Color', 'b')
        xlabel('TIME [s]')
        ylabel('r_' + string(j) + ' [m]')
    end
    sgtitle('Position. Inertial frame')
    
    
    figure()
    for j = 1:3
        subplot(3, 1, j)
        plot(time, state_true(:, j+3), 'LineWidth', 2, 'Color', 'g')
        xlabel('TIME [s]')
        ylabel('v_' + string(j) + ' [m/s]')
    end
    sgtitle('Velocity. Inertial frame')

    figure()
    for j = 1:3
        subplot(3, 1, j)
        plot(time, stateDot_true(:, j+3), 'LineWidth', 2, 'Color', ...
            "#FF00FF")
        xlabel('TIME [s]')
        ylabel('a_' + string(j) + '[m/s^2]')
    end
    sgtitle('Acceleration. Inertial frame')
end

