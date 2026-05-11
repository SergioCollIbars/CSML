clear;
clc;
close all;

load("/Users/sergiocollibars/Desktop/DLA_Gravity_Gradiometry/Simulator/data/output.mat");
load("/Users/sergiocollibars/Desktop/DLA_Gravity_Gradiometry/Simulator/data/time.mat");
load("/Users/sergiocollibars/Desktop/DLA_Gravity_Gradiometry/Simulator/data/NB_MOON.mat");

data = mC_struct.sim1;
Pf = data.Pf;
Xf = data.Xf;

BN_mat = compute_orientation_SC(time, Xf(1:6,:)', "RTN", NB_MOON_mat);

C_pos_bias_x = nan(length(Pf(:, 1)), 12);
C_pos_bias_y = nan(length(Pf(:, 1)), 12);
C_pos_bias_z = nan(length(Pf(:, 1)), 12);

C_bias_SF   = nan(length(Pf(:, 1)), 6, 6);

sigma_pos    = nan(length(Pf(:, 1)), 3);
for k = 1:length(Pf(:, 1))
    p = reshape(Pf(k, :), [18, 18]);

    maxInd = 3*k;
    minInd = maxInd - 2;
    BN = BN_mat(minInd:maxInd, :);

    A   = blkdiag(BN, BN, eye(6),eye(6));
    P_B = A * p * A';

    sigma = sqrt(diag(P_B));        % standard deviations
    C = P_B ./ (sigma * sigma');    % correlation matrix

    C_pos_bias_x(k, :) = C(1, 7:end);
    C_pos_bias_y(k, :) = C(2, 7:end);
    C_pos_bias_z(k, :) = C(3, 7:end);

    C_bias_SF(k, 1, :) = C(7, 13:end);
    C_bias_SF(k, 2, :) = C(8, 13:end);
    C_bias_SF(k, 3, :) = C(9, 13:end);
    C_bias_SF(k, 4, :) = C(10, 13:end);
    C_bias_SF(k, 5, :) = C(11, 13:end);
    C_bias_SF(k, 6, :) = C(12, 13:end);

    sigma_pos(k, :)    = sigma(1:3);
end

figure()
plot(1:length(Pf(:, 1)), C_pos_bias_x, 'LineWidth', 2); grid on;
xlabel('Time [s]'); title('X pos direction');
legend('b_{xx}', 'b_{xy}', 'b_{xz}', 'b_{yy}', 'b_{yz}', 'b_{zz}', ...
    'SF_{xx}', 'SF_{xy}', 'SF_{xz}', 'SF_{yy}', 'SF_{yz}', 'SF_{zz}');


figure()
plot(1:length(Pf(:, 1)), C_pos_bias_y, 'LineWidth', 2); grid on;
xlabel('Time [s]'); title('Y pos direction');
legend('b_{xx}', 'b_{xy}', 'b_{xz}', 'b_{yy}', 'b_{yz}', 'b_{zz}', ...
    'SF_{xx}', 'SF_{xy}', 'SF_{xz}', 'SF_{yy}', 'SF_{yz}', 'SF_{zz}');

figure()
plot(1:length(Pf(:, 1)), C_pos_bias_z, 'LineWidth', 2); grid on;
xlabel('Time [s]'); title('Z pos direction');
legend('b_{xx}', 'b_{xy}', 'b_{xz}', 'b_{yy}', 'b_{yz}', 'b_{zz}', ...
    'SF_{xx}', 'SF_{xy}', 'SF_{xz}', 'SF_{yy}', 'SF_{yz}', 'SF_{zz}');

figure()
plot(1:length(Pf(:, 1)), 3.*sigma_pos, 'LineWidth', 2); grid on;
xlabel('Time [s]'); title('3 \sigma position');
legend('radial', 'tangetial', 'normal');

for k = 1:6
    figure();
    plot(1:length(Pf(:, 1)), squeeze(C_bias_SF(:, k, :)), 'LineWidth', 2); 
end
grid on;
xlabel('Time [s]'); title('bias correlation to SF');
legend('SF_{xx}', 'SF_{xy}', 'SF_{xz}', 'SF_{yy}', 'SF_{yz}', 'SF_{zz}');