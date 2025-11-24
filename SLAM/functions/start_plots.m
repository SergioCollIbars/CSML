function [fig, h, ax] = start_plots(time, n_max, Xc_t, ...
    Xc_0, P0c_grav, Nbatch)
%=== EKF Summary Figure: All Plots in One Window ===%

% Assumes variables:
%   time   - time vector in seconds
%   Ncs    - number of gravity coefficients
% You can adapt as needed.

t_hours = time(:).'/3600;   % time in hours for x-axis limits

idx = 1:Nbatch:length(t_hours);
xv  = t_hours(idx);
grayColor = [.7 .7 .7];

fig = figure('Name','EKF Summary: State, Bias, Gravity',...
    'NumberTitle','off');

tiledlayout(fig, 3, 3, 'TileSpacing','compact', 'Padding','compact');

% Preallocate handle structs
h   = struct();  % line handles
ax  = struct();  % axes handles

%% 1) Position error (top-left)
nexttile; hold on; grid on;
h.Err_pos = plot(nan, nan, '-', 'LineWidth', 1.5);  % error curve
h.Cov_pos = plot(nan, nan, '-', 'LineWidth', 1.5);  % 3σ curve
xlabel('Time [hr]');
ylabel('[m]');
title('Position Error Norm');
xlim([0, t_hours(end)]);
ax.pos = gca;
for k = 2:length(xv)
    xline(xv(k), '-', 'LineWidth', 0.3, 'Color', 'k');
end

%% 2) Velocity error (top-middle)
nexttile; hold on; grid on;
h.Err_vel = semilogy(nan, nan, '-', 'LineWidth', 1.5);
h.Cov_vel = semilogy(nan, nan, '-', 'LineWidth', 1.5);
xlabel('Time [hr]');
ylabel('[mm/s]');
title('Velocity Error Norm');
xlim([0, t_hours(end)]);
ax.vel = gca;
set(ax.vel, 'YScale','log');
for k = 2:length(xv)
    xline(xv(k), '-', 'LineWidth', 0.3, 'Color', 'k');
end

%% 3–8) Bias plots (middle + bottom, 6 tiles)
h.BiasErr = gobjects(1, 6);
h.BiasCov = gobjects(1, 6);
ax.bias   = gobjects(1, 6);

for k = 1:6
    nexttile; hold on; grid on;
    h.BiasErr(k) = plot(nan, nan, '-', 'LineWidth', 1.5);
    h.BiasCov(k) = plot(nan, nan, '-', 'LineWidth', 1.5);
    title(sprintf('Bias %d', k));
    xlabel('Time [hr]');
    ylabel('[mE]');
    xlim([0, t_hours(end)]);
    ax.bias(k) = gca;
    for j = 2:length(xv)
        xline(xv(j), '-', 'LineWidth', 0.3, 'Color', 'k');
    end
end

%% 9) Gravity field coefficients (bottom-right)
nexttile; hold on; grid on;
h.Err_c = semilogy(nan, nan, '-', 'LineWidth', 1.5); % |error|
h.Cov_c = semilogy(nan, nan, '-', 'LineWidth', 1.5); % 3σ
xlabel('Coeff index');
ylabel('[-]');
title('Gravity Coefficient Error');
xlim([2, n_max+1]); ylim([1E-7, 1E-1]);
xticks(1:n_max);
ax.grav = gca;
set(ax.grav, 'YScale','log');

% compute initial gravity error
errC   = abs(Xc_0 - Xc_t);
sigmaC = sqrt(diag(P0c_grav));
[X_true, xvals] = orderValues(abs(Xc_t(2:end)), n_max);
[errC_order, ~] = orderValues(errC(2:end), n_max);
[sigmaC_order, ~] = orderValues(sigmaC(2:end), n_max);

semilogy(xvals, X_true, 'Color', 'k', 'LineWidth', 2)
set(h.Err_c, 'XData', xvals, 'YData',...
    errC_order, 'Color', 'r', 'Marker','*',...
    "MarkerSize", 2, 'LineStyle','none');
set(h.Cov_c, 'XData', xvals, 'YData', 3.*sigmaC_order, ...
    'LineWidth', 2, 'Color', 'r');

%=== At this point you have:
% fig  -> figure handle
% h    -> struct of line handles (Err_pos, Cov_pos, Err_vel, Cov_vel,
%           BiasErr(1:6), BiasCov(1:6), Err_c, Cov_c)
% ax   -> struct of axes handles (pos, vel, bias(1:6), grav)
%
% You can update any plot later with, for example:
%   set(h.Err_pos, 'XData', t_hours, 'YData', pos_err_norm);
%   set(h.Cov_pos, 'XData', t_hours, 'YData', 3*sigma_pos_norm);
%   set(h.BiasErr(3), 'XData', t_hours, 'YData', bias3_err);
%   set(h.Err_c, 'XData', 1:Ncs-1, 'YData', abs(errC(2:end)));
%   set(h.Cov_c, 'XData', 1:Ncs-1, 'YData', 3*sigmaC(2:end));

end
