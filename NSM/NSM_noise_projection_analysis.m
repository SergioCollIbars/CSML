clear;
clc;
close all;
format long g;

addpath('functions/')
addpath('../QGG_gravEstim/src/')
addpath('../QGG_navigation/data/')
set(0,'defaultAxesFontSize',16);

%%      EVALUATE NOISE PROJECTION
% Description: Undertand and evaluate noise projection in Null space
% method.
% Author: Sergio Coll
% Date: 05/21/25

% select error type
error = "position"; % options: position / attitude

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]

poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.32E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 1;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
PHI0 = reshape(eye(6,6), [36, 1]);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC, 0), t, [r0;v0; PHI0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';
plot_trajectory(state_t, "BENNU");

% compute noise sensitivity
Hmat  = ones(length(t), 6) * NaN;
% % [Hmat] = compute_Noise_Sensitivity(t, rn, vn, RA, DEC, W0, W, Cnm, Snm, Re, GM, n_max, normalized, error);

%  compute residual projecction
[coeffs_vals, terms] = compute_coeff_proj(t, rn, vn, RA, DEC, W0, W, Cnm, Snm, Re, GM, n_max, normalized, error);

% compute information loss after the projection
[I_proj, I_org, I_proj_ideal] = compute_Inf_loss(t, rn, vn, RA, DEC, W0, W, Cnm, Snm, Re, GM, n_max, normalized, error);

% Plot results
[num_C, num_S, str_C, str_S] = SH_xlabel(n_max);

% Account for 00 term
num_C = num_C + 1;
num_C = [1, num_C];
str_C = ["C_{00}", str_C];

figure()
B = 1 - diag(I_proj./I_org); % information ratio
subplot(1, 2, 1)
plot(1:Nc, B(1:Nc), 'LineWidth', 2)
xticks(num_C);
xticklabels(str_C);
grid on;
title('C_{nm}')
subplot(1, 2, 2)
plot(1:Ns, B(Nc+1:Ncs), 'LineWidth', 2)
xticks(num_S);
xticklabels(str_S);
grid on;
title('S_{nm}')
sgtitle('Information loss due to projection')

figure()
B = diag(I_proj./I_proj_ideal);           % projected information ratio
subplot(1, 2, 1)
semilogy(1:Nc, B(1:Nc), 'LineWidth', 2)
xticks(num_C);
xticklabels(str_C);
grid on;
title('C_{nm}')
subplot(1, 2, 2)
semilogy(1:Ns, B(Nc+1:Ncs), 'LineWidth', 2)
xticks(num_S);
xticklabels(str_S);
grid on;
title('S_{nm}')
sgtitle('Projected information ratio w.r.t spherical case')

for j = 1:length(coeffs_vals(1, 1, :))
    figure()
    plot(t, coeffs_vals(:, :, j)', 'LineWidth', 2)
    if(length(terms) == 6)
        legend('dY_{RR}','dY_{RT}', 'dY_{RN}', 'dY_{TT}', 'dY_{TN}', 'dY_{NN}')
    else
        legend('dY_{RT}', 'dY_{RN}', 'dY_{TT}', 'dY_{TN}', 'dY_{NN}')
    end
    
    title('$dY \cdot \vec{v}_{' + string(j) + '}$', 'Interpreter','latex');
end

for j = 1:3
    figure()
% %     y = Hmat(:, :, j);
    y = smoothdata(Hmat(:, :, j), 1, 'movmean', 10);
    plot(t, real(y), 'LineWidth', 2)
    if(error == "position")
        legend('', 'S_{\sigma^2_{RT}}', 'S_{\sigma^2_{RN}}',...
            'S_{\sigma^2_{TT}}', 'S_{\sigma^2_{TN}}', 'S_{\sigma^2_{NN}}');
    else
        legend('S_{\sigma^2_{RR}}', 'S_{\sigma^2_{RT}}', 'S_{\sigma^2_{RN}}',...
            'S_{\sigma^2_{TT}}', 'S_{\sigma^2_{TN}}', 'S_{\sigma^2_{NN}}');
    end
end
% % legend('','S_{\sigma^2_{xy}}', 'S_{\sigma^2_{xz}}',...
% %     'S_{\sigma^2_{yy}}', 'S_{\sigma^2_{yz}}', 'S_{\sigma^2_{zz}}');


%% FUNCTIONS
function [coeffs_vals, terms] = compute_coeff_proj(t, rn, vn, RA, DEC, W0, W, Cnm, Snm, Re, GM, n_max, normalized, error)
    syms y11 y12 y13 y22 y23 y33 real
    coeffs_vals = ones(length(t), 6, 3) * NaN;
    for k = 1:length(t)
        disp('Running .. ' + string(k/length(t)));
        
        % position and velcoity
        r  = rn(:, k); v = vn(:, k);
    
        % rotation to RTN frame
        ACI_RTN = RTN2ECI(r, v);
    
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % rotation to body frame
        ACAF_BODY = ACAF_ACI * ACI_RTN;
% %         ACAF_RTN = ACAF_ACI;
    
        % position / attitude partials & null space
        if(error == "position")
            [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
                 ACAF_ACI, ACAF_BODY);
             C = null([Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
             dY = [y12 y13 y22 y23 y33]';

% %              C = null([Hp(1,:);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
% %              dY = [y11 y12 y13 y22 y23 y33]';
        elseif(error == "attitude")
            [Hp] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_BODY);
            C = null([Hp(1,:); Hp(2,:);Hp(3,:);Hp(5, :); Hp(6, :); Hp(9, :)]');
            dY = [y11 y12 y13 y22 y23 y33]';
        else
            warning('Select a correct mode for the error type')
            break;
        end
        Nterms  = length(dY);
        Nspace  = Nterms - 3;

        % meas residual proyection
        dY_proj  = C' * dY;
    
        % get the coefficients
        for j = 1:Nspace
            [coeffs_vals(k, 1:Nterms, j), terms] = coeffs(dY_proj(j), dY);
        end
    end
end

function [Hmat] = compute_Noise_Sensitivity(t, rn, vn, RA, DEC, W0, W, Cnm, Snm, Re, GM, n_max, normalized, error)
    % noise values
    Hmat = ones(length(t), 6 , 3) * NaN;
    
    % create noise values
    syms v11 v12 v13 v21 v22 v23 v31 v32 v33 positive
    R = diag([v11,v12,v13,v22,v23,v33]);

    for k = 1:length(t)
        disp('Running .. ' + string(k/length(t)));
        
        % position and velcoity
        r  = rn(:, k); v = vn(:, k);
    
        % rotation to RTN frame
        ACI_RTN = RTN2ECI(r, v);
    
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % rotation to body frame
        ACAF_RTN = ACAF_ACI * ACI_RTN;
% %         ACAF_RTN = ACAF_ACI;
        if(error == "position")
            [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
            ACAF_ACI, ACAF_RTN); 
            C = null([Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
            R = diag([v12,v13,v22,v23,v33]);
        elseif(error == "attitude")
            [Hp] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_RTN);
            C = null([Hp(1, :); Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :); Hp(9, :)]');
            R = diag([v11, v12,v13,v22,v23,v33]);
        else
            warning('Select a correct error mode...');
        end
        % noise projection
        r = C' * R * C;
    
        [~, D] = eig(r);
        d = diag(D);
        
        for i = 1:length(d)
            % compute the sensitivity matrix
            J11 = diff(d(i), v11);
            J12 = diff(d(i), v12);
            J13 = diff(d(i), v13);
            J22 = diff(d(i), v22);
            J23 = diff(d(i), v23);
            J33 = diff(d(i), v33);
        
            % evaluate around a nominal sensitivity value
            val = 1E-10;
            H11 = subs(J11, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H12 = subs(J12, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H13 = subs(J13, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H22 = subs(J22, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H23 = subs(J23, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
            H33 = subs(J33, [v11, v12, v13, v22, v23, v33], [val, val, val, val, val, val]);
        
            Hmat(k, :, i)  = [H11, H12, H13, H22, H23, H33];
        end
    end
end

function [I_proj, I_org, I_proj_ideal] = compute_Inf_loss(t, rn, vn, RA, DEC, W0, W, Cnm, Snm, Re, GM, n_max, normalized, error)
    I_proj = 0; I_org = 0; I_proj_ideal = 0;
    sigma_RR = 1; sigma_RT = 1; sigma_RN = 1;
    sigma_TT = 1; sigma_TN = 1; sigma_NN = 1;
    for k = 1:length(t)
        disp('Running .. ' + string(k/length(t)));
        
        % position and velcoity
        r  = rn(:, k); v = vn(:, k);
    
        % rotation to RTN frame
        ACI_RTN = RTN2ECI(r, v);
    
        % ACAF to ACI rotation matrix
        Wt = W0 + W * t(k);
        ACAF_ACI =rotationMatrix(pi/2 + RA, pi/2 - DEC, Wt, [3, 1, 3]);
    
        % rotation to body frame
        ACAF_BODY = ACAF_ACI * ACI_RTN;
        ACI_BODY  = ACI_RTN;
% %         ACAF_BODY = ACAF_ACI;
% %         ACI_BODY = eye(3,3);

        % compute gravity field partials
        planetParams = [GM, Re, n_max, normalized];
        [~, ~, H] = gradiometer_meas(t(k) ,planetParams, ACAF_ACI, [r', v'], ...
                zeros(9, length(t)), Cnm, Snm, ACI_BODY);
    
        % position partials & null space
        if(error == "position")
            [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
                 ACAF_ACI, ACAF_BODY);
             C = null([Hp(2,:);Hp(3,:);Hp(5,:);Hp(6, :);Hp(9, :)]');
             R = diag([sigma_RT, sigma_RN, sigma_TT,sigma_TN,sigma_NN].^2);
             R_id = eye(5,5);
             h = [H(4,:);H(7,:);H(5,:);H(8, :);H(9,:)];
             
% %              C  = null([Hp(3,:);Hp(5,:);Hp(6,:);Hp(9, :)]');
% %              R = diag([sigma_RN,sigma_TT, sigma_TN,sigma_NN].^2);
% %              R_id = eye(4,4);
% %              h = [H(7,:);H(5,:);H(8,:);H(9, :)];

             h_org = [H(1, :);H(4,:);H(7,:);H(5, :);H(8, :);H(9, :)];
             R_org = diag([sigma_RR,sigma_RT, sigma_RN, sigma_TT, sigma_TN, sigma_NN].^2);
% %              h_org = h;
% %              R_org = R;
             
        elseif(error == "attitude")
            [Hp] = compute_rotPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ACAF_ACI, ACAF_BODY);
            C = null([Hp(1,:); Hp(2,:);Hp(3,:);Hp(5, :); Hp(6, :); Hp(9, :)]');
            R = diag([sigma_RR, sigma_RT, sigma_RN, sigma_TT, sigma_TN, sigma_NN].^2);
            R_id = eye(6,6);
            h = [H(1, :);H(4,:);H(7,:);H(5, :);H(8, :); H(9, :)];


% %             C = null([Hp(1,:);Hp(5,:);Hp(6, :);Hp(9,:)]');
% %             R = diag([sigma_RR,sigma_TT,sigma_TN, sigma_NN].^2);
% %             R_id = eye(4,4);
% %             h = [H(1,:);H(5,:);H(8,:);H(9, :)];
% % 
            h_org = [H(1, :);H(4,:);H(7,:);H(5, :);H(8, :); H(9, :)];
            R_org = diag([sigma_RR,sigma_RT, sigma_RN, sigma_TT, sigma_TN, sigma_NN].^2);
% %             h_org = h;
% %             R_org = R;
        else
            warning('Select a correct mode for the error type')
            break;
        end

        % project matrices
        h_proj = C' * h;
        r = C' * R * C;
        r_ideal = C' * R_id * C;

        % compute information
        I_org = I_org + h_org' * inv(R_org) * h_org;
        I_proj = I_proj + h_proj' * inv(r) * h_proj;
        I_proj_ideal = I_proj_ideal + h_proj' * inv(r_ideal) * h_proj;
    end
end