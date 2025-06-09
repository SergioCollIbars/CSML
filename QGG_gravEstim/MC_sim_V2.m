clear;
clc;
close all;
set(0,'defaultAxesFontSize',16);

%%              MONTE CARLO SIMULATION
% Inputs
MC = 10;         % number of MC runs

% outputs
GM1 = zeros(MC, 2);
B1  = zeros(MC, 6);
D1  = zeros(MC, 6); 

% position & velocity errors
errP = ones(1, 3) * 0;      % [m]
errV = ones(1, 3) * 0;        % [m/s]

% run MC simulation
for mc = 1:MC
    disp("MC simulation: " + string(mc) + "/" + string(MC))
    % run main
    run("main.m")

    % read acc and orbit data
    T = readtable('accData.txt');
    T2 = readtable('orbitData.txt');
    dataAcc = table2array(readtable('accData.txt'));
    dataOrbit = table2array(readtable('orbitData.txt'));

    TIME = dataAcc(:, 1);
    N = length(TIME);
    AT = TIME(2) - TIME(1);
    G = dataAcc(:, 2:10);
    r_ACI = dataOrbit(:, 2:4);
    v_ACI = dataOrbit(:, 8:10);

    % include position and velocity errors
    % % r_ACI = r_ACI + ones(length(r_ACI(:, 1)), 3).*errP;
    r_ACI = r_ACI + normrnd(0, 0.07, [length(r_ACI(:, 1)), 3]);
    v_ACI = v_ACI + ones(length(v_ACI(:, 1)), 3).*errV;

    sdiag  = 10*6.32E-12;
    sndiag = 10*2.52E-10;
    % % sndiag = sdiag;
    R = diag([sdiag, sndiag, sndiag, sdiag, sndiag, sdiag].^2);

    % Do estimation
    [X1, ID, P, ~] = estimateGM(N, AT, G, r_ACI, R);
    GM1(mc, 1) = X1(1);
    GM1(mc, 2) = 5.2 - X1(1);
    B1(mc, :) = abs([X1(2), X1(4), X1(6), X1(8), X1(10), X1(12)]).*1E-10;
    D1(mc, :) = abs([X1(3), X1(5), X1(7), X1(9), X1(11), X1(13)]).*1E-10;
        
    clc;
end

% % % estimate R
% % [R_new] = update_R_val(R, N, AT, G, r_ACI);

% compute statistics
m = mean(GM1(:, 1));
s = std(GM1(:, 1));

% plot estimation value.
figure()
errorbar(1, m, s, 'LineWidth', 2, 'Marker', 'sq')
xlim([-0.5, 3.5])
ylim([5.19, 5.2])
title('GM estimation value')
set(gca, 'XTick',1, 'XTickLabel',"Actual")
grid on;

% compute statistics
B = [-282.1677,8.9215,16.4064,-636.5541,-19509.3451,-113.2609] * 1E-9 * 2;
D = [0.0065,-0.0103,-0.0001,0.2029,0.9381,0.0033] * 1E-9/(86400)*2;
mB = mean(B1);
sB = std(B1);
mD = mean(D1);
sD = std(D1);

% formal uncertainty
S = diag(sqrt(P));
sB_F = [S(2); S(4); S(6); S(8); S(10); S(12)].*1E-10;
sD_F = [S(3); S(5); S(7); S(9); S(11); S(13)].*1E-10;

% dispal SNR in micro-Gal (50cm acc baseline)
SNR_B = B./sB;
SNR_D = D./sD;
disp('SNR bias  :' + string(SNR_B'));
disp('SNR drift :' + string(SNR_D'));

SNR_B_F = B'./sB_F;
SNR_D_F = D'./sD_F;
disp('Formal SNR bias  :' + string(SNR_B_F));
disp('Formal SNR drift :' + string(SNR_D_F));

% plot estimation value.
figure()
semilogy(1:6, B, 'k-', 'LineWidth', 1.5)
hold all;
errorbar(1:6, mB, sB, 'LineWidth', 2, 'Color', 'b', 'Marker', 'sq')
title('Bias estimation value')
legend('Truth', 'estimate 1 \sigma')
grid on;

figure()
semilogy(1:6, D, 'k-', 'LineWidth', 1.5)
hold all;
errorbar(1:6, mD, sD, 'LineWidth', 2, 'Color', 'b', 'Marker', 'sq')
title('Drift estimation value')
set(gca,'YScale','log')
legend('Truth', 'estimate 1 \sigma')
grid on;

figure()
semilogy(TIME./86400, abs(ID), 'LineWidth', 2);
xlabel('TIME [days]')
legend('Inertial', 'RTN')
title('Information matrix determinant')

%%                      FUNCTIONS
function [X1,ID, P, r] = estimateGM(N, AT, G, r_ACI, R)
    ID = ones(2, N) * NaN;
    std_theta = 10 * pi/(180*3600);    
    A_theta   = normrnd(0, std_theta,  [3, N]);      % [rad]

    % generate measurements. Gamma RTN
    G_RTN = ones(3*N, 3)*NaN;
    for k = 1:N
        % rotation matrix
        ACI_RTN = eye(3, 3);
        ACI_BODY = rotationMatrix(A_theta(1, k), A_theta(2, k), A_theta(3, k), [3, 2, 1]);
        up = 3*k;
        down = up -2;
    

        % GM effect (ACI)
        r3 = vecnorm(r_ACI(k, :))^3;
        r5 = vecnorm(r_ACI(k, :))^5;
        x = r_ACI(k, 1);
        y = r_ACI(k, 2);
        z = r_ACI(k, 3); 
        H_GM = 5.2.*[1/r3-3*x*x/r5; -3*x*y/r5; -3*x*z/r5; 1/r3-3*y*y/r5;...
            -3*z*y/r5; 1/r3-3*z*z/r5];
        G_GM = [H_GM(1), H_GM(2), H_GM(3);...
            H_GM(2), H_GM(4), H_GM(5);...
            H_GM(3), H_GM(5), H_GM(6)];

        % Gamma in ACI frame
        G_ACI = [G(k, 1), G(k, 2), G(k, 3);...
            G(k, 4), G(k, 5), G(k, 6);...
            G(k, 7), G(k, 8), G(k, 9)];

        % Remove GM effect
        G_B = G_ACI - G_GM;

        % include GM with attitutde error
        G_GM_B = ACI_BODY' * G_GM * ACI_BODY;
        G_ACI = G_B + G_GM_B;
        G_RTN(down:up, :) =  ACI_RTN' * G_ACI * ACI_RTN;
    end
    
    % LS estimate biases and drifts. Using Gamma
    Ax = 0; Ax2  = 0;
    Nx = 0;
    Rinv = inv(R);
    for k = 1:N
        r3 = vecnorm(r_ACI(k, :))^3;
        r5 = vecnorm(r_ACI(k, :))^5;
        x = r_ACI(k, 1);
        y = r_ACI(k, 2);
        z = r_ACI(k, 3); 

        up =3*k;
        down = up -2;
        g = G_RTN(down:up, :);
        yy = [g(1,1);g(1,2);g(1,3);g(2,2);g(2,3);g(3,3)];
        i = k - 1;
        
        h_GM = [1/r3-3*x*x/r5; -3*x*y/r5; -3*x*z/r5; 1/r3-3*y*y/r5;...
            -3*z*y/r5; 1/r3-3*z*z/r5];
        h2_GM = 1/r3*[-2; 0; 0; 1; 0; 1];
        h_b  = zeros(6, 2*6);
        for j = 1:6
            maxPos = 2*j;
            minPos = maxPos - 1;
            h_b(j, minPos:maxPos) = [1E-10, 1E-10 * i *AT];
        end
        h  = [h_GM, h_b];
        h2 = [h2_GM, h_b];

    
        % increment Normal equation and information matrix
        Ax  = Ax + h' * Rinv * h;       % cartesian coord
        Ax2 = Ax2 + h2' * Rinv * h2;    % spherical coord
        Nx  = Nx + h' * Rinv * yy;
% % 
% %         % Information Determinant (ID)
% %         p1 = inv(Ax); 
% %         p2 = inv(Ax2);
% %         ID(1, k) = det(p1(1, 1));
% %         ID(2, k) = det(p2(1, 1));
    end
    
    % solve state
    X1 = Ax\Nx;
    P = inv(Ax);

    % compute residuals
    r = ones(6, N) * NaN;
    for k = 1:N
        up =3*k;
        down = up -2;
        g = G_RTN(down:up, :);
        yy = [g(1,1);g(1,2);g(1,3);g(2,2);g(2,3);g(3,3)];
        i = k - 1;
        
        h_GM = [1/r3-3*x*x/r5; -3*x*y/r5; -3*x*z/r5; 1/r3-3*y*y/r5;...
            -3*z*y/r5; 1/r3-3*z*z/r5];
        h_b  = zeros(6, 2*6);
        for j = 1:6
            maxPos = 2*j;
            minPos = maxPos - 1;
            h_b(j, minPos:maxPos) = [1E-10, 1E-10 * i *AT];
        end
        h  = [h_GM, h_b];
        
        r(:, k) = yy - h * X1;
    end
end


function [R] = update_R_val(R_init, N, AT, G, r_ACI)
% Initialization
R = R_init;  % Initial guess for R (6 x 6)
max_iter = 10;
tol = 1e-26;
alpha = 0.8;

for k = 1:max_iter
    % Run batch filter with current R to get residuals (6 x Nt)
    [~,~, ~, residuals] = estimateGM(N, AT, G, r_ACI, R);
    
    % Compute sample covariance of residuals
    Nt = size(residuals, 2);
    r_mean = mean(residuals, 2);
    residuals_centered = residuals - r_mean;
    R_sample = (residuals_centered * residuals_centered') / (Nt - 1);
    
    % Update R
    R_new = alpha * R_sample + (1 - alpha) * R;
    
    % Check convergence
    if norm(R_new - R, 'fro') < tol
        disp(['Converged at iteration ', num2str(k)]);
        break;
    end
    
    % Update
    R = R_new;
end

% Now R is the re-estimated measurement noise covariance matrix
end