clear;
clc;
close all;

% Description: script to create plot showing the error asociated to the
% orbit radius, noise measurement and delta GM.

% Input
GM = 5.2;
r = linspace(245, 3000, 1000);

% Matrices
sigma_AR_noise = zeros(1, length(r));
sigma_AR_GM = zeros(1, length(r));

% loop
for k = 1:length(r)

    rt = r(k);

    % compute error
    sigma_AR_noise(k) = (rt^4) / (3 * GM);
    sigma_AR_GM(k) = rt / (3 * GM);
end


% plot
figure();

semilogy(r, sigma_AR_noise, 'LineWidth', 1.5, 'Color', 'k');
hold on;
semilogy(r, sigma_AR_GM, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle','--');
grid on;
set(gca, 'FontSize', 15);
xlabel('Orbit radius [m]', 'FontSize', 12);
ylabel('\sigma (\Delta_R/\epsilon_i) [m]', 'FontSize', 12);

title('Uncertanty in the radial postion created by model errors', 'FontSize', 14);
legend('Noise effect', '\delta GM effect');

