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
%GM = 5.2;
J2 = 1.08E-3;                               % J2 radious

Re = 6563E3;                                % Earth radious

Nm = 9;
Ns = 1;

Nf = 10;
Nv = 100;
Nr = Nf;

% Initial conditions. Orbital elements
e  = 0;
a  = Re + 250E3;
rho = a.*(1 - e^2);
f = 0;

i = deg2rad(96);
Omega = 0;
omega = 0;

% time values
n = sqrt(GM./(a.^3));
T = 10 * pi./n;

t_max = T;
t_min = 0;

% sigma values
frec = linspace(2E-4, 0.1, Nf);                 % frec range [Hz]
sigmaM = linspace(1E-4, 1E-3, Nv);              % sensor variance [E  * Hz-0.5]
%sigmaM = linspace(1E-3, 1, Nv);                % sensor variance [E  * Hz-0.5]
meanValue = 0;

% variables definiton
err = zeros(Nf, Nv);
Xvar = zeros(Nf, Nv);

for j =1:Nf
    for s = 1:Nv
        disp("frec number: " + string(j) + "/" + Nf);
        disp("  variance number: " + string(s) + "/" + Nv);

        % current values
        rho_k = rho;
        a_k = a;
        frec_k = frec(j);
        tmax_k = t_max;
        sigmaM_k = sigmaM(1);

        % compute initial state vectors (orbit frame)
        r0 = rho_k / (1 + e * cos(f)) * [cos(f);...
                                  sin(f);...
                                  0];
        v0 = sqrt(GM / rho_k) * [-sin(f);...
                          e + cos(f);...
                          0];
    
        % compute rotation matrix
        R = rotationMatrix(Omega, i, omega, [3,1,3]);

        % rotate initial state vectors. {I, J, K} frame
        r0 = R' * r0;
        v0 = R' * v0;

        % compute number of points
        Nt = round(tmax_k * frec_k);
        t = linspace(t_min, tmax_k, Nt);

        % solve orbit
        X0 = [r0(1); r0(2); r0(3); v0(1); v0(2); v0(3)];
        
        options = odeset('RelTol',1e-13,'AbsTol',1e-13);
        [~, state] = ode113(@(t, x) EoM(t, x, GM, J2, Re, Ns), t, X0, options);
        state = state';

        % compute sigma value
        sigma = sigmaM_k * sqrt(frec_k) * 1E-9;       % [s^-2]

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
        
        % time loop
        for k = 1:Nt      
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

        % error
        if(Ns == 2)
            Xtrue = [GM, J2]';
        else 
            Xtrue = GM; 
        end
        err(j, s) = abs(X - Xtrue);
        Xvar(j, s) = sqrt(inv(Ax));
        
        clc;
    end
end

% GOCE GM error
GMerr = 1E-5 * a(1)^3;

% NASA GM error
GMerr2 = 1E-15 * a(1)^3;

% plot error per variance 
figure();
for k =1:Nf
% %     plot(sigmaM, err(k, :), 'o' , 'DisplayName', "frec " + string(frec(k)) + ...
% %         " Hz");
    loglog(sigmaM, err(k, :), 'o' );
    hold all;
end
%plot(sigmaM, GMerr * ones(1, Nv), '--', 'DisplayName', "GOCE GM error req");
plot(sigmaM, GMerr2 * ones(1, Nv), '--', 'DisplayName', "NASA GM error req");

grid on;
xlabel("noise variance [E / Hz^{0.5}]");
ylabel("X err [m^3 / (s^2)]");
title("Estimation error per noise variance. Monte Carlo run N = " + string(Nf));


% plot error per frecuency 
figure();
for k =1:Nv
% %     plot(sigmaM, err(k, :), 'o' , 'DisplayName', "frec " + string(frec(k)) + ...
% %         " Hz");
    loglog(frec, err(:, k), 'o' );
    hold all;
end
%plot(sigmaM, GMerr * ones(1, Nv), '--', 'DisplayName', "GOCE GM error req");
plot(frec, GMerr2 * ones(1, Nf), '--', 'DisplayName', "NASA GM error req");

grid on;
xlabel("measurement frequency [Hz]");
ylabel("X err [m^3 / (s^2)]");
title("Estimation error per frecuency. Monte Carlo run N = " + string(Nv));
legend("\sigma_v = " + string(sigmaM(1)) + "E/Hz^{0.5}")


% plot state variance per variance 
figure();
for k =1:Nf
    plot(sigmaM, Xvar(k, :),'DisplayName', "frec " + string(frec(k)) + ...
        " Hz");
    hold all;
end
grid on;
xlabel("noise variance [E / Hz^{0.5}]");
ylabel("estate variance [1/s^2]");
title("Estimation variance per noise at different frequencies");