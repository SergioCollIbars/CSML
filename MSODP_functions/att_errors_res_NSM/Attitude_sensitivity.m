clear;
clc;
close all;
addpath('functions/');
set(0,'defaultAxesFontSize',16);

%%                 ATTITUDE SENSITIVITY ANALYSIS
% Description: Evaluate the attitude sensitivity for the GG measurements
% per gradiometer component.
% Author: Sergio Coll-Ibars
% Date: 06-24-2026
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Extract GG observations
folderPath = "/Users/sergiocollibars/Documents/GG_observations/120by120/500km/";
[GG_obs]   = parser_GG_obs_MSODP(folderPath); Nt = length(GG_obs(:, 1));
t0         = datetime(2000,1,1,12,0,0,'TimeZone','UTC');
t_dateTime = t0 + seconds(GG_obs(:, 1));

%% Compute sensitivity
disp('Computing attitude sensitivity ...');
sensitivity = nan(6, 3, Nt);

att_Err           = zeros(3, Nt);
[GG_rot, GG_nom]  = rotate_GG_obs(GG_obs,att_Err);
mask              = [1, 1, 1, 0, 1, 1, 0, 0, 1];
for j = 1:Nt
    [Hrot_1st]           = compute_rotPartials_analy(GG_nom(j, 2:end)', ...
                            eye(3));
    sensitivity(:, :, j) = Hrot_1st(logical(mask), :);   

end

%% plot results
error_att_title = ["yaw", "pitch", "roll"];
component_title = ["xx", "xy", "xz", "yy", "yz", "zz"];
t_initial       = datetime(2008,8,1,0,0,0, 'TimeZone','UTC');
t_final         = datetime(2008,8,2,0,0,0, 'TimeZone','UTC');
arcsecond2rad   = pi / (180 * 3600); % [rad / arcseconds]
for comp = 1:6
    figure();
    for i = 1:3
        subplot(1, 3, i);
        data = squeeze(sensitivity(comp, i, :)).* (1/1E-15*arcsecond2rad); 
        % [micro-E/arcsecond]
        plot(t_dateTime, data, 'LineWidth', 1.2, ...
            'Color', 'k');
        title(error_att_title(i)); grid on;
        xlim([t_initial t_final]); 
        ylabel('\mu E / arcsecond');
    end
    sgtitle("Gradiometer component " + component_title(comp));
end