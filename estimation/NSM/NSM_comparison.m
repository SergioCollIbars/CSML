clear;
clc;
close all;
format long g;
addpath('../functions/')
addpath('../../../QGG_gravEstim/src/')
set(0,'defaultAxesFontSize',16);

%%              NSM METHODS COMPARISON
% Description: Compare Rummel's formulation vs the null space approach.
% Test in small radius orbits using Linear Leas Squares.
% Author: Sergio Coll
% Date: 09/28/24

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

% % path = "HARMCOEFS_EROS_CD_1.txt";
% % [Cnm, Snm, Re] = readCoeff(path);
% % n_max  = 10;
% % normalized = 1;
% % GM =  459604.431484721;          % Point mass value    [m^3/s^2]
% % W = 1639.38928 * pi/180 /86400;  % Rotation ang. vel   [rad/s]
% % W0 = 0;                          % Initial asteroid longitude
% % RA = deg2rad(11.363);            % Right Ascension     [rad]
% % DEC = deg2rad(17.232);           % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.3E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% position error
Ar = 5*[1;1;1];            % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 3;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% noise values
noise0 = zeros(9, Nt);
sigma  = 6.32E-12;
noise  =  normrnd(0, sigma, [9, length(t)]);

% Integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
STM0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0;STM0], options);
rn = state_t(:, 1:3)' + ones(3, Nt).*Ar;
vn = state_t(:, 4:6)';

% generate measurements
[Y, ~, ~] = gradiometer_meas(t ,asterParams, poleParams, state_t, ...
                noise0, Cnm, Snm);

% perturb nominal coefficient
[X] = mat2list(Cnm, Snm, Nc, Ns);
sigma_n = [1E-5;1E-5;1E-5;1E-6;1E-6];
[Xp] = perturb_coeff(sigma_n, n_max, X);
[Cp, Sp] = list2mat(n_max, Nc, Ns, Xp);

% Gravity estimation
Ax_R = 0; Ax_N = 0;
Nx_R = 0; Nx_N = 0;
R_R = diag([sqrt(2)*sigma, sigma].^2);
R_N = diag([sigma, sigma, sigma, sigma, sigma].^2);
n = ones(2, Nt)*NaN;
n_trunc = n_max; 
for j = 1:Nt
    % computed meas.
    [~, Hc_ACI, Hc_RTN] = gradiometer_meas(t(j) ,asterParams, poleParams, [rn(:, j)', vn(:, j)'], ...
            noise0, Cnm, Snm);
    
    % RTN rotation matrix
    ACI_RTN = RTN2ECI(rn(:, j), vn(:, j));
    rn_RTN = ACI_RTN' * rn(:, j);
    rn_ACI = rn(:, j);
    
    % ACAF to ACI rotation matrix
    Wt = W0 + W * t(j);
    ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);

     % compute Point mass approx error
    [Hpos] = compute_posPartials(n_trunc, normalized, Cp, Sp, Re, GM, rn_RTN, ACAF_ACI*ACI_RTN);

    % Rummel's method
    [ax, nx] = rummels_method(Y(:, j), Hc_RTN, R_R, ACI_RTN, noise(:, j));
    Ax_R  = Ax_R + ax;
    Nx_R  = Nx_R + nx;

    % Null space method
    [ax, nx] = nullSpace_method(Y(:, j), Hc_RTN, R_N, Hpos, ACI_RTN, noise(:, j));
    Ax_N  = Ax_N + ax;
    Nx_N  = Nx_N + nx;
end

% solve LS
SH_R = Ax_R\Nx_R;
sh_N = Ax_N\Nx_N;
SH_N = sh_N(2:end);

P_N =  inv(Ax_N);
P_R =  inv(Ax_R);
sigma_N = sqrt(diag(P_N));
sigma_R = sqrt(diag(P_R));


[num_C, num_S, str_C, str_S] = SH_xlabel(n_max);
[nc, ns, ~] = count_num_coeff(n_trunc); 
[bias_C, bias_S] = compute_bias(n_trunc);

% plot trajectory
plot_trajectory(state_t, 'BENNU');
figure()
plot(t./T, vecnorm(state_t(:, 1:3)'), 'LineWidth', 2)
hold on;
plot(t./T, ones(1, Nt)*Re, 'LineWidth', 2, 'Color', 'r', 'LineStyle','--')
xlabel('Orb. Period number, T')
ylabel('[m]')
title('Orbit radius along trajectory. T = ' + string(T./3600) + ' h')
legend('orbit radius', 'brillouin sphere')

% plot SH estimation
figure()
subplot(1, 2, 1)
semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold all;
semilogy(1:Nc-1, abs(SH_R(1:Nc-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
semilogy(1:Nc-1, abs(SH_N(1:Nc-1)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
title('Estimation SH Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;

subplot(1, 2, 2)
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold on;
semilogy(1:Ns, abs(SH_R(Nc:Nc+Ns-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
semilogy(1:Ns, abs(SH_N(Nc:Nc+Ns-1)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','auto')
title('Estimation SH. Snm coefficients')
xticks(num_S);
xticklabels(str_S);
grid on;
legend('truth','n=0 truncation', 'null-space')

% plot uncertainty
figure()
subplot(1, 2, 1)
semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold all;
semilogy(1:Nc-1, abs(sigma_R(1:Nc-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'm', 'MarkerFaceColor', 'auto')
semilogy(1:Nc-1, abs(sigma_N(1:Nc-1)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'g', 'MarkerFaceColor', 'auto')
title('Uncertainty Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;

subplot(1, 2, 2)
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold on;
semilogy(1:Ns, abs(sigma_R(Nc:Nc+Ns-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'm', 'MarkerFaceColor', 'auto')
semilogy(1:Ns, abs(sigma_N(Nc:Nc+Ns-1)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'g', 'MarkerFaceColor','auto')
title('Uncertainty Snm coefficients')
xticks(num_S);
xticklabels(str_S);
grid on;
legend('truth','n=0 truncation', 'null-space')

% ploting difference
figure()
subplot(1, 2, 1)
semilogy(1:Nc-1, abs(X(2:Nc)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold all;
semilogy(1:Nc-1, abs(X(2:Nc) - SH_R(1:Nc-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
semilogy(1:Nc-1, abs(X(2:Nc) - SH_N(1:Nc-1)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor', 'auto')
title('Estimation difference SH Cnm coefficients')
xticks(num_C);
xticklabels(str_C);
grid on;


subplot(1, 2, 2)
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'k', 'MarkerFaceColor', 'auto')
hold on;
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns) - SH_R(Nc:Nc+Ns-1)), 'Marker','square', 'LineStyle','-', 'LineWidth', 2, 'Color', 'r', 'MarkerFaceColor', 'auto')
semilogy(1:Ns, abs(X(Nc+1:Nc+Ns) - SH_N(Nc:Nc+Ns-1)), 'Marker','diamond', 'LineStyle','-', 'LineWidth', 2, 'Color', 'b', 'MarkerFaceColor','auto')
title('Estimation SH. Snm coefficients')
xticks(num_S);
xticklabels(str_S);
grid on;
legend('truth','n=0 truncation', 'null-space')
%% FUNCTIONS
function [ax, nx] = rummels_method(Y, Hc, R, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1), Y(2), Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;
    
    % select measurements
    dY = [dy(5) - dy(9); dy(8)];
    hc = [Hc(5, 2:end) - Hc(9, 2:end); Hc(8, 2:end)];
    
    % information and normal matrices
    ax = hc' * inv(R) * hc;
    nx = hc' * inv(R) * dY;
end

function [ax, nx] = nullSpace_method(Y, Hc, R, Hp, ACI_RTN, noise)
    % reshape meas [ACI]
    ddU = [Y(1),Y(2),Y(3);Y(4),Y(5),Y(6);Y(7),Y(8),Y(9)];
 
    % Rotate meas to RTN coordinates
    dy_RTN = ACI_RTN' * ddU * ACI_RTN;
    dy = reshape(dy_RTN, [9, 1]) + noise;

    % select measurements
    dY = [dy(1);dy(4);dy(7);dy(5);dy(8)];

    % look for null space
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');

    % project measurements
    y  = C' * dY;
    hc = C' * [Hc(1, 1:end); Hc(4, 1:end); Hc(7, 1:end);Hc(5, 1:end);...
        Hc(8, 1:end)];
    r  = C' * R * C;

    % information and normal matrices
    ax = hc' * inv(r) * hc;
    nx = hc' * inv(r) * y;
end

function [num_C, num_S, str_C, str_S] = SH_xlabel(n_max)    
    num_C = ones(1, n_max-1) * NaN;
    num_S = num_C;

    str_C = cell(1, n_max - 1);
    str_S = str_C;

    num_C(1) = 1;
    for j = 3:n_max
        num_C(j-1) = j + num_C(j-2);
    end
    
    num_S(1) = 1;
    for j = 3:n_max
        num_S(j-1) = (j-1) + num_S(j-2);
    end

    for j = 2:n_max
        str_C{j - 1} = "C_{" + string(j) + "0}";
        str_S{j - 1} = "S_{" + string(j) + string(j) + "}";
    end 
end

function [Xp] = perturb_coeff(sigma_n, n_max, X)
    [Nc, Ns, ~] = count_num_coeff(n_max); 

    % perturbed values
    Xp = X;

    m  = 0;
    n = 2;
    for j =2:Nc
        Xp(j) = Xp(j) + normrnd(0, sigma_n(n-1));
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 0;
        end
    end
    m  = 1;
    n = 2;
    for j =Nc+1:Nc+Ns
        Xp(j) = Xp(j) + normrnd(0, sigma_n(n-1));
        if(m < n)
            m = m + 1;
        else
            n = n + 1;
            m = 1;
        end
    end
end

function [bias_C, bias_S] = compute_bias(n_max)
    % SH harmonics
    [Nc, Ns, ~] = count_num_coeff(n_max); 
    
    % bias vector (GM not considered)
    bias_C = ones(Nc-1, 1) * NaN;
    bias_S = ones(Ns, 1)   * NaN;
    
    % fill n = 2 values
    bias_C(1:3) = 2/3;
    bias_S(1:2) = 2/3;
    
    % fill C bias values
    val = 4;
    for j = 3:n_max
        maxInd = j + val;
        minInd = val;
        bias_C(minInd:maxInd) = j/3;

        val = maxInd + 1;
    end
    
    % fill S bias values
    val = 3;
    for j = 3:n_max
        maxInd = j - 1 + val;
        minInd = val;
        bias_S(minInd:maxInd) = j/3;

        val = maxInd + 1;
    end
end