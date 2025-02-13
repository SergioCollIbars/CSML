clear;
close all;
clc;

format long g;

%%           PERFECT QGG DEVICE COEFFICENT ESTIMATION

% Decription: computes the GM coefficient for a 2BP assuming a perfect QGG
%   device plus some noise.

% Imports
addpath('functions/');

% frame
frame = "inertial";                         % body / inertial

% Inputs
GM = 3.9833E14;                             % Gravity param [m^3 s^-2]
J2 = 1.08E-3;

Re = 6563E3;                                % Earth radious

r0 = [7563E3, 0, 0]';                       % Initial position vector [m];

v0 = [0, 7257.28966892064, 0]';  % i = 0 degree
v0 = [0 ,5131.67873792886, 5131.67873792886]';  % i = 45 degree
v0 = [0 ,0, 7257.28966892064]';  % i = 90 degree
%%v0 = [500, 7257.28966892064, 0]'; % i = 0 degree eliptical orbit

Nt = 1000000;
Nm = 9;
Ns = 2;

n = sqrt(GM/  (r0(1)^3));
T = 10 * pi / n;

Nt = 100000;

t_max = T;
t_min = 0;

t = linspace(t_min, t_max, Nt);

sigma = 1E-6;                                  % Sensor variance [m]
meanValue = 0;

Nf = 10;
Nv = 100;

frec = linspace(2E-2, 1E-1, Nf);                % frec range Hz
sigmaM = linspace(8E-3, 11E-3, Nv);             % sensor variance E Hz-1

% solve orbit
X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];

options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state] = ode113(@(t, x) EoM(t, x, GM, J2, Re, Ns), t, X0, options);
state = state';

% plot 3D orbit
figure();
plot3(state(1, :), state(2, :), state(3, :));
grid on;
%axis equal;



% variables definition
Y = zeros(Nm, Nt);
Hi = zeros(Nm, Nt*Ns);

dP = zeros(Nm, Nt);

% Initial normal equation values
Ax = 0;
Nx = 0;

% values to save
detBN = zeros(1, Nt);
U = zeros(3, Nt);


for k = 1:Nt
    %disp(" Iteration: " + string(k) + "/" + string(Nt));
    
    % Inertial coords
    xi = state(1, k);
    yi = state(2, k);
    zi = state(3, k);
    
    % compute potential
    u = potentialGradient_GM(GM, xi, yi, zi);
    if(Ns == 2)
        u = u + potentialGradient_J2(GM, J2, Re, xi, yi, zi);
    end
    
    % get observables
    if(Ns == 2)
        accT = -potentialGradient2_GM(GM, xi, yi, zi) + ...
            -potentialGradient2_J2(GM, J2, Re, xi, yi, zi);
    else
        accT = - potentialGradient2_GM(GM, xi, yi, zi);
    end
    
    % express in body or inertial
    if frame == "body"
        ri = state(1:3, k);
        vi = state(4:6, k);
    
        alpha = orbitalElem(ri, vi, GM);
        
        omega = alpha(8);
        Omega = alpha(7);
        i = alpha(6);
        a = alpha(3);
    
        if(isnan(Omega))
            R1 = rotationMatrix(0, i, 0, [3, 1, 3]);
        else
            R1 = rotationMatrix(Omega, i, omega, [3, 1, 3]);
        end
        
        ro = R1 * ri;
        vo = R1 * vi;
        
        h = cross(ro, vo);
    
        omega_k = cross(ro, vo) / (vecnorm(ro)^2);
        omega_k_u = omega_k / vecnorm(omega_k);
    
        Omega_k = -2 * vecnorm(h) / (vecnorm(ro)^4) * ...
                dot(ro, vo) * omega_k_u;
        Uo = R1 * u;
        %%Omega_k = cross(ro, Uo);
        
        theta_k = atan2(ro(2), ro(1));
    
        R2 = rotationMatrix(theta_k, 0, 0, [3, 1, 3]);
    
        % rotation matrix
        BN = R2 * R1;
        
        
        % compute square ang velocity
        Momega2 = [-(omega_k(2)^2 + omega_k(3)^2), omega_k(1) * omega_k(2), omega_k(1) * omega_k(3);...
                   omega_k(1) * omega_k(2), -(omega_k(1)^2 + omega_k(3)^2), omega_k(2) * omega_k(3);...
                   omega_k(1) * omega_k(3), omega_k(2) * omega_k(3), -(omega_k(1)^2 + omega_k(2)^2)];
    
        % compute angular acc matrix
        MOmega = [0, -Omega_k(3), Omega_k(2);...
                  Omega_k(3), 0, -Omega_k(1);...
                  -Omega_k(2), Omega_k(1), 0];
    
        % add angular values to measurements. Body frame
        accT = BN * accT * BN' + MOmega + Momega2;
    else
        Momega2 = zeros(3, 3);
        MOmega = zeros(3, 3);
    
        BN = eye(3, 3);
    end
    
    % save potential gradient
    U(:, k) = BN * u;
    
    % save determinant BN
    detBN(k) = det(BN);
    
    % add noise to data
    accT = accT + sigma * randn(size(accT)) + meanValue;
    
    Ydata = reshape(accT, [Nm, 1]) - reshape(MOmega, [Nm, 1]) - ...
        reshape(Momega2, [Nm, 1]);
    Y(:, k) = Ydata;
    
    % compute measurement sensitivity matrix
    HGM_GM = -BN * potentialGradient2_GM(GM, xi, yi, zi) * BN';
    HGM_GM = reshape(HGM_GM./GM, [Nm, 1]);
    
    HJ2_GM = -BN * potentialGradient2_J2(GM, J2, Re, xi, yi, zi) * BN';
    HJ2_GM = reshape(HJ2_GM./GM, [Nm, 1]);
    
    HGM_J2 = zeros(Nm, 1);
    
    HJ2_J2 = -BN * potentialGradient2_J2(GM, J2, Re, xi, yi, zi) * BN';
    HJ2_J2 = reshape(HJ2_J2./J2, [Nm, 1]);
    
    HJ2_B = -BN * potentialGradient2_J2(GM, J2, Re, xi, yi, zi) * BN';
    HJ2_B = reshape(HJ2_B./(J2*GM), [Nm, 1]);
    
    if(Ns ==  2)
        % % H = [(HGM_GM + HJ2_GM),(HGM_J2 + HJ2_J2)];
        H = [HGM_GM, HJ2_B];
    else
        H = HGM_GM;
    end
    
    % stack visibility matrix
    up = k * Ns;
    down = up - (Ns - 1);
    Hi(:, down:up) = H;
    
    
    % run Batch filter
    sigma2_R = sigma^2;
    Ri = eye(Nm, Nm) * sigma2_R;
    
    [Ax, Nx] = batchFilter(Ydata, H, Ri, Ax, Nx);
        
end
        
% solve normal eq.
X = Ax\ Nx;
if(Ns == 2), X(2) = X(2)/X(1); end

disp("Estimation value [km^3 s^2]");
x = X;
x(1) = x(1) * 1E-9;
disp(x);


% error
disp("Error in km^3 s^2: ")
if(Ns == 2)
    Xtrue = [GM, J2]';
else 
    Xtrue = GM; 
end
%err(1) = X - Xtrue;

% % err(1) = err(1) * 1E-9;
% % disp(err);

Ytilde = [];
sigmaX0 = 0;

% compute postfit
for k = 1:Nt
    up = k * Ns;
    down = up - (Ns - 1);

    y = Y(:, k);
    hi = Hi(:, down:up);
    dp =  y - hi * X;
    dP(:, k) = reshape(dp, [Nm, 1]);

    Ytilde = [Ytilde; y];
    sigmaX0 = sigmaX0 + hi' * inv(Ri) * hi;
end

sigmaX0 = sqrt(inv(Ax));
disp("STD of the estimate: " + string(sigmaX0 * 1E-9));

% variance of the data
sigmaY = std(Ytilde);
disp("STD of the data: " + string(sigmaY));

% error sample variance
sigmaSerr = sigmaY / sqrt(Nt);
disp(" Error sample variance: " + string(sigmaSerr));

% plot postfit
figure();

for k = 1:Nm
    subplot(3, 3, k);
    
    plot(t, dP(k, :), '*');
    MdP = mean(dP(k, :));

    xlabel("TIME [s]");
    ylabel("dP_" + string(k));
    title("mean: " + string(MdP));
    grid on;
end

% plot determinant BN
figure();

plot(t, detBN);
xlabel("TIME");
ylabel("|BN|");
title("Rotation matrix determinant");
grid on;

% plot gravity potential gradient
figure();

subplot(3, 1, 1);
plot(t, U(1, :));
xlabel("TIME");
ylabel("\Delta U_x");
grid on;

subplot(3, 1, 2);
plot(t, U(2, :));
xlabel("TIME");
ylabel("\Delta U_y");
grid on;

subplot(3, 1, 3);
plot(t, U(3, :));
xlabel("TIME");
ylabel("\Delta U_z");

sgtitle("Potential gradient in " + frame + " coordinates");
grid on;



