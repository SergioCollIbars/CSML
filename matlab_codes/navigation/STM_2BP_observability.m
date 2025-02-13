clear;
clc;
close all;
format long g;

addpath('functions/')
set(0,'defaultAxesFontSize',16);
%%          2BP STM PROPAGATOR
% Description: Test the relative measurements observability using the 
% STM propagation method.

% Simulation conditions
G = 6.67430e-11;    % [N m^2 Kg^-2]
N = 40;

% define planet parameters
R =  6.3781E6;  % [m]
m = 5.9722E24;  % [Kg]
GM = G*m;       % [m^3 s^-2]
n = sqrt(GM / R^3); % mean motion circular orbit [1/s]
nmax = 3;
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

frec = 1/10;   % [Hz]
Nt =  2*pi/n*frec;
TIME = linspace(0, 2*pi/n, Nt);
% TIME = [0, 2*pi/n];

% Initital conditions. Orbital parameters
e = 0.4;                     % eccentrycity
a = R + 500e3;               % semi major axis [m]

i = deg2rad(0);             % inclination [rad]
omega = deg2rad(0);          % arg periapsis [rad]
Omega = deg2rad(0);          % RAAN [rad]
f = deg2rad(-180);           % true anomaly [rad]

e_array = linspace(0.1, 0.8, N);        % [-]
a_array = linspace(R, R + 20E4, N);     % [m]
i_array = linspace(0, pi/2, N);         % [rad]

[Xe, Xa] = meshgrid(e_array, a_array);

% absolute state uncertainty
p_abs = ones(N, N) * NaN;
v_abs = ones(N, N) * NaN;

p_rel = ones(N, N) * NaN;
v_rel = ones(N, N) * NaN;
for j = 1:N
    disp('j = ' + string(j))
    for h = 1:N
        e = Xe(h, j);
        a = Xa(h, j);
        rho = a * (1 - e^2);   % orbital param [m]
        
        % solve initial cond
        [r0, v0] = orbElems_2_ACI(rho, f, GM, Omega, omega, i, e);
        
        % define initial conditions. Small body
        X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
        
        % integatre trajectory and STM
        STM_0 = reshape(eye(6,6), [36, 1]);
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        [t, state_true] = ode113(@(t, x) EOM_navigation(t, x, planetParams, ...
            poleParams, Cmat, Smat, "2BP"), TIME, [X0; STM_0], options);

        % compute uncertainty
        O = zeros(6*length(t), 12);
        STM = state_true(:, 7:end);
        H = [eye(3,3), zeros(3,3);zeros(3,3), eye(3,3)];
        for  k = 1:length(t)-1
            PHI_a = reshape(STM(k, :), 6,6);
            PHI_t = reshape(STM(k+1, :), 6,6);
            maxInd = 6*k;
            minInd  = maxInd - 5;
            O(minInd:maxInd, :) = [H * (PHI_t - PHI_a) , 0.5*(H* (PHI_t + PHI_a))];
            %O(minInd:maxInd, :) = [];
        end
        if(rank(O) == 12)
            l  = O'*O;
            A  = l(1:6, 1:6);
            B  = l(1:6, 7:12);
            C  = l(7:12, 1:6);
            D  = l(7:12, 7:12);
            
            P11 = inv(A - B*inv(D)*C);
            P12 = -P11*B*inv(D);
            P21 = -inv(D)*C*P11;
            P22 = inv(D) + inv(D)*C*P11*B*inv(D);
            
            P = [P11, P12;P21, P22];
            sigma = diag(P);
            
            p_abs(h, j) = vecnorm(sigma(1:3));
            v_abs(h, j) = vecnorm(sigma(4:6));
            p_rel(h, j) = vecnorm(sigma(7:9));
            v_rel(h, j) = vecnorm(sigma(10:12));
        else
            w = 1;
        end
    end
end

figure()
subplot(2, 2, 1)
contourf(Xe, Xa./1000, p_abs)
xlabel('e [-]')
ylabel('a [km]')
colorbar()
title('Abs pos \sigma')
set(gca,'ColorScale','log');
subplot(2, 2, 2)
contourf(Xe, Xa./1000, v_abs)
xlabel('e [-]')
ylabel('a [km]')
colorbar()
title('Abs vel \sigma')
set(gca,'ColorScale','log');

subplot(2, 2, 3)
contourf(Xe, Xa./1000, p_rel)
xlabel('e [-]')
ylabel('a [km]')
colorbar()
title('Rel pos \sigma')
set(gca,'ColorScale','log');
subplot(2, 2, 4)
contourf(Xe, Xa./1000, v_rel)
xlabel('e [-]')
ylabel('a [km]')
colorbar()
title('Rel vel \sigma')
set(gca,'ColorScale','log');




