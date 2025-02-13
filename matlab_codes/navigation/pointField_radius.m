clc;
clear;
close all;
format long g;

%%          RESOLVE RADIUS IN POINT FIELD


simulation = "CR3BP_inertial"; % options: 2BP / CR3BP

% define planet parameters
G = 6.67430e-11;    % [N m^2 Kg^-2]
if(simulation == "2BP")
    R =  6.3781E6;  % [m]
    m = 5.9722E24;  % [Kg]
    GM = G*m;       % [m^3 s^-2]
    n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
    planetParams = [GM, R, 0, 0];
    poleParams = [0, 0, 0, 0];

    % Simulation conditions
    tmin = 0;           % [s]
    tmax = 1/n;         % [s]
    frec = 10/10;       % [Hz]
    N = tmax*frec;
    t = linspace(tmin, tmax, N);
    nmax = 0;

elseif(simulation == "CR3BP_inertial")
    m_1 = 5.974E24;  % [Kg]
    m_2 = 7.348E22;  % [Kg]
    GM = G*(m_1 + m_2); % [m^3 s^-2]
    planetParams = [m_2/(m_1+m_2); 0; 0; 0];

    % Simulation conditions
    tmin = 0;           % [s]
    tmax = 100;         % [s]
    frec = 10/10;       % [Hz]
    N = tmax*frec;
    t = linspace(tmin, tmax, N);
    nmax = 0;
end


if(simulation == "2BP")
    % plot potential and eigen values direction
    plot_potential_eigen_2BP(planetParams, poleParams, -R, R);
elseif(simulation == "CR3BP")
    % plot potential and eigen values direction
    plot_potential_eigen_3CRBP(planetParams, -1.5, 1.5)
elseif(simulation == "CR3BP_inertial")
    % plot gradiometer measurements
    plot_grad_CR3BP(planetParams, -1, 1)
end

if(simulation == "2BP")
    % Initital conditions. Orbital parameters
    e = 0.4;                     % eccentrycity
    a = R + 100e3;               % semi major axis [m]
    rho = a * (1 - e^2);         % orbital param [m]
    
    i = deg2rad(45);             % inclination [rad]
    omega = deg2rad(0);          % arg periapsis [rad]
    Omega = deg2rad(0);          % RAAN [rad]
    f = deg2rad(-180);            % true anomaly [rad]
    
    [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
    
    % define initial conditions. Small body
    X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
end

% % % integatre trajectory and STM
% % STM_0 = reshape(eye(6,6), [36, 1]);
% % 
% % % define integration options
% % options = odeset('RelTol',1e-13,'AbsTol',1e-13);
% % 
% % % ODE 113
% % [time, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
% %     poleParams, Cmat, Smat, "2BP"), t, [X0; STM_0], options);
% % 
% % % reconstruct trajectory
% % gamma = zeros(9, N);
% % for j =1:N
% %     [~, ddU_ACI] = gradiometer_meas(t(j), state_true(j,1:3)', ...
% %     Cmat, Smat, poleParams, planetParams);
% %     gamma(:, j) = reshape(ddU_ACI, [9, 1]);
% % end



%% FUNCTIONS

function plot_potential_eigen_2BP(planetParams, poleParams, rmin, rmax)
    % create meshgrid
    N = 25;
    x = linspace(rmin, rmax, N);
    y = linspace(rmin, rmax, N);
    z = 0;
    [X, Y] = meshgrid(x, y);
    GM = planetParams(1);
    Cmat = [1, 0, 0, 0;...
     0, 0, 0, 0;...
    -1.08E-3, 0, 1.57E-6, 0;...
     2.53E-6, 2.18E-6, 3.11E-7, 1.02E-7]; 

    % potential function
    U = zeros(N, N);
    v1x = zeros(N, N);
    v1y = zeros(N, N);
    for j = 1:N % y direction
        for i = 1:N % x direction
            rvec = [X(i, j); Y(i, j); z];
            r = vecnorm(rvec);
            U(i, j) = - GM/r;
            [~, T] = gradiometer_meas(1, rvec, ...
                Cmat, zeros(3,3), poleParams, planetParams);
            T  =reshape(T, [3,3]);
            [V, D] = eig(T);
            [~, p] = max(max(D));    % max eigenvalue

            v1x(i, j) = V(1, p);
            v1y(i, j) = V(2, p);
        end
    end

    % plot potential
    figure()
    pcolor(X, Y, U)
    xlabel('X [m]')
    ylabel('Y [m]')
    zlabel('U')
    zlim([-10, 0])
    title('2BP gravity potential')
    shading interp
    hold on;
    quiver(X, Y, v1x, v1y, 'Color', 'r', 'LineWidth', 2)
end


function plot_grad_CR3BP(planetParams, rmin, rmax) % plot_potential_eigen_3CRBP
    % create meshgrid
    N = 50;
    x = linspace(rmin, rmax, N);
    y = linspace(rmin, rmax, N);
    z = 0;
    [X, Y] = meshgrid(x, y);
    mu = planetParams(1);

    eps = 0.01;

    % potential function
    U11 = ones(N, N)*NaN;
    U12 = U11;
    U13 = U11;
    U22 = U11;
    U23 = U11;
    U33 = U11;
    for j = 1:N % y direction
        for i = 1:N % x direction
            rvec = [X(i, j); Y(i, j); z];
            r1 = vecnorm(rvec - [-mu;0;0]);
            r2 = vecnorm(rvec - [1-mu;0;0]);

            [ddU_I] = gradmeas_CR3BP_inertial(mu, rvec(1), ...
                rvec(2), rvec(3), 0);

            % save values
            if((r1 > eps) && (r2>eps))
                U11(i, j) = ddU_I(1,1);
                U12(i, j) = ddU_I(1,2);
                U13(i, j) = ddU_I(1,3);
                U22(i, j) = ddU_I(2,2);
                U23(i, j) = ddU_I(2,3);
                U33(i, j) = ddU_I(3,3);
            end
        end
    end

    % plot potential
    figure()
    subplot(2, 3, 1)
    surf(X, Y, U11)
    shading interp
    ylabel('U_{11}')

    subplot(2, 3, 2)
    surf(X, Y, U12)
    shading interp
    ylabel('U_{12}')

    subplot(2, 3, 3)
    surf(X, Y, U13)
    shading interp
    ylabel('U_{13}')

    sgtitle('CR3BP gradiometer values. t=0')
end

function plot_potential_eigen_3CRBP(planetParams, rmin, rmax)
    % create meshgrid
    N = 50;
    x = linspace(rmin, rmax, N);
    y = linspace(rmin, rmax, N);
    z = 0;
    [X, Y] = meshgrid(x, y);
    mu = planetParams(1);

    % potential function
    U = zeros(N, N);
    v1x = zeros(N, N);
    v1y = zeros(N, N);
    for j = 1:N % y direction
        for i = 1:N % x direction
            rvec = [X(i, j); Y(i, j); z];
            r1 = rvec - [-mu;0;0];
            r2 = rvec - [1-mu;0;0];
            U(i, j) = - mu./vecnorm(r2) - (1-mu)./vecnorm(r1);
            T = gradmeas_rotFrame(mu,rvec(1), rvec(2), rvec(3), ...
                vecnorm(r1), vecnorm(r2));
            T = reshape(T, [3,3]);
            [V, D] = eig(T);
            [~, p] = max(max(abs(D)));    % max eigenvalue

            phi = atan2(V(3, 3), V(3, 2));
            lambda = atan2(V(2, 3), V(1, 3));

            v1x(i, j) = cos(lambda);
            v1y(i, j) = sin(lambda);
        end
    end

    % plot potential
    figure()
    surf(X, Y, U)
    xlabel('X [-]')
    ylabel('Y [-]')
    zlabel('U')
    zlim([-10, 0])
    title('CR3BP gravity potential')
    shading interp
    hold on;
    quiver(X, Y, v1x, v1y, 'Color', 'r', 'LineWidth', 2)
end


