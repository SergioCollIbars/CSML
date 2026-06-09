% https://isdc-data.gfz.de/grace/Level-2/GFZ/RL06/
%% read_all_GSM_compute_mean_and_C20_trend.m
clear; clc; close all;

%% User inputs
folder_path = '/Users/sergiocollibars/Documents/GRACE_GSM_files/';

n_max = 96;   % Use 96 for the file you showed
R_E = 6378136.3; %#ok<NASGU> % Earth reference radius [m], optional

%% Find GSM files
files = dir(fullfile(folder_path, 'GSM-*'));

Nfiles = length(files);

if Nfiles == 0
    error('No GSM files found in folder: %s', folder_path);
end

%% Preallocate coefficient arrays
% Dimensions: degree x order x time
C_all      = nan(n_max+1, n_max+1, Nfiles);
S_all      = nan(n_max+1, n_max+1, Nfiles);
sigmaC_all = nan(n_max+1, n_max+1, Nfiles);
sigmaS_all = nan(n_max+1, n_max+1, Nfiles);

time_decimal = nan(Nfiles,1);
time_datetime = NaT(Nfiles,1);

%% Read all files
for k = 1:Nfiles

    filename = fullfile(folder_path, files(k).name);

    fprintf('Reading file %d/%d: %s\n', k, Nfiles, files(k).name);

    [Cnm, Snm, sigmaC, sigmaS] = read_GRACE_GSM_file(filename, n_max);

    C_all(:,:,k)      = Cnm;
    S_all(:,:,k)      = Snm;
    sigmaC_all(:,:,k) = sigmaC;
    sigmaS_all(:,:,k) = sigmaS;

    % Extract approximate time from filename, e.g.
    % GSM-2_2015102-2015131_GRAC_GFZOP_BB01_0600
    [time_datetime(k), time_decimal(k)] = extract_GSM_time(files(k).name);

end

%% Sort files by time
[time_decimal, idx_sort] = sort(time_decimal);
time_datetime = time_datetime(idx_sort);

C_all      = C_all(:,:,idx_sort);
S_all      = S_all(:,:,idx_sort);
sigmaC_all = sigmaC_all(:,:,idx_sort);
sigmaS_all = sigmaS_all(:,:,idx_sort);

%% Compute coefficient means
C_mean      = mean(C_all, 3, 'omitnan');
S_mean      = mean(S_all, 3, 'omitnan');
sigmaC_mean = mean(sigmaC_all, 3, 'omitnan');
sigmaS_mean = mean(sigmaS_all, 3, 'omitnan');

%% Save mean matrices
save(fullfile(folder_path, 'GSM_coefficient_means.mat'), ...
    'C_mean', 'S_mean', 'sigmaC_mean', 'sigmaS_mean', ...
    'time_datetime', 'time_decimal', 'n_max');

fprintf('\nSaved mean coefficients to GSM_coefficient_means.mat\n');

%% Plot C20 time evolution and trend
n = 2;
m = 0;

C20 = squeeze(C_all(n+1,m+1,:));

valid = ~isnan(C20) & ~isnat(time_datetime);

t = time_decimal(valid);
C20_valid = C20(valid);
time_valid = time_datetime(valid);

% Linear trend: C20(t) = a*t + b
p = polyfit(t, C20_valid, 1);

C20_trend = polyval(p, t);

trend_per_year = p(1);

fprintf('\nC20 linear trend = %.6e per year\n', trend_per_year);

figure;
plot(time_valid, C20_valid, 'o-', 'LineWidth', 1.2);
hold on;
plot(time_valid, C20_trend, 'k--', 'LineWidth', 2);
grid on;

xlabel('Time');
ylabel('$\bar{C}_{20}$', 'Interpreter', 'latex');
title('$\bar{C}_{20}$ Time Evolution and Linear Trend', 'Interpreter', 'latex');

legend( ...
    '$\bar{C}_{20}$ monthly values', ...
    sprintf('Linear trend = %.3e / year', trend_per_year), ...
    'Interpreter', 'latex', ...
    'Location', 'best');

%% Optional: plot C20 anomaly relative to the mean
C20_mean = C_mean(n+1,m+1);
C20_anomaly = C20_valid - C20_mean;

p_anom = polyfit(t, C20_anomaly, 1);
C20_anom_trend = polyval(p_anom, t);

figure;
plot(time_valid, C20_anomaly, 'o-', 'LineWidth', 1.2);
hold on;
plot(time_valid, C20_anom_trend, 'k--', 'LineWidth', 2);
grid on;

xlabel('Time');
ylabel('$\Delta \bar{C}_{20}$', 'Interpreter', 'latex');
title('$\bar{C}_{20}$ Anomaly Relative to Mean', 'Interpreter', 'latex');

legend( ...
    '$\Delta \bar{C}_{20}$', ...
    sprintf('Linear trend = %.3e / year', p_anom(1)), ...
    'Interpreter', 'latex', ...
    'Location', 'best');

%% Compute Hydrology + Ice truth curve from GRACE GSM anomalies

R_E = 6378136.3;  % Earth reference radius [m]

Nt = size(C_all, 3);

hydro_ice_geoid_mm = nan(n_max+1, Nt);

for k = 1:Nt

    % Monthly anomaly relative to the mean field
    dC = C_all(:,:,k) - C_mean;
    dS = S_all(:,:,k) - S_mean;

    for n = 0:n_max

        degree_sum = 0;

        for m = 0:n

            degree_sum = degree_sum + ...
                dC(n+1,m+1)^2 + dS(n+1,m+1)^2;

        end

        hydro_ice_geoid_mm(n+1,k) = 1000 * R_E * sqrt(degree_sum);

    end

end

%% RMS Hydrology + Ice curve over all months

hydro_ice_rms_mm = sqrt(mean(hydro_ice_geoid_mm.^2, 2, 'omitnan'));

degrees = 0:n_max;

figure;
semilogy(degrees, hydro_ice_rms_mm, 'k:', 'LineWidth', 2);
grid on;

xlabel('Spherical harmonic degree');
ylabel('Geoid signal [mm]');
title('Hydrology + Ice Truth Signal from GRACE GSM');
legend('Hydrology + Ice RMS truth', 'Location', 'best');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local function: read one GSM file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Cnm, Snm, sigmaC, sigmaS] = read_GRACE_GSM_file(filename, n_max)

    Cnm    = nan(n_max+1, n_max+1);
    Snm    = nan(n_max+1, n_max+1);
    sigmaC = nan(n_max+1, n_max+1);
    sigmaS = nan(n_max+1, n_max+1);

    fid = fopen(filename, 'r');

    if fid < 0
        error('Could not open file: %s', filename);
    end

    while ~feof(fid)

        line = fgetl(fid);

        if ~ischar(line)
            continue;
        end

        line_trim = strtrim(line);

        if startsWith(line_trim, 'GRCOF2')

            data = textscan(line_trim, '%s %d %d %f %f %f %f %s %s %s');

            n = data{2};
            m = data{3};

            if n <= n_max && m <= n_max

                Cnm(n+1,m+1)    = data{4};
                Snm(n+1,m+1)    = data{5};
                sigmaC(n+1,m+1) = data{6};
                sigmaS(n+1,m+1) = data{7};

            end

        end

    end

    fclose(fid);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local function: extract time from GSM filename
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [date_mid, decimal_year] = extract_GSM_time(fname)

    % Example filename:
    % GSM-2_2015102-2015131_GRAC_GFZOP_BB01_0600
    %
    % This extracts:
    % start date = year 2015, day-of-year 102
    % end date   = year 2015, day-of-year 131
    % and uses the midpoint as the time tag.

    expr = '_(\d{4})(\d{3})-(\d{4})(\d{3})_';
    tokens = regexp(fname, expr, 'tokens');

    if isempty(tokens)
        warning('Could not extract date from filename: %s', fname);
        date_mid = NaT;
        decimal_year = nan;
        return;
    end

    tokens = tokens{1};

    year_start = str2double(tokens{1});
    doy_start  = str2double(tokens{2});
    year_end   = str2double(tokens{3});
    doy_end    = str2double(tokens{4});

    date_start = datetime(year_start,1,1) + days(doy_start - 1);
    date_end   = datetime(year_end,1,1)   + days(doy_end - 1);

    date_mid = date_start + (date_end - date_start)/2;

    year_mid = year(date_mid);

    date_year_start = datetime(year_mid,1,1);
    date_year_end   = datetime(year_mid+1,1,1);

    decimal_year = year_mid + days(date_mid - date_year_start) / ...
        days(date_year_end - date_year_start);

end