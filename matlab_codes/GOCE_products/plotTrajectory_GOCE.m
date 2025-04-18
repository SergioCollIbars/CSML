clear;
clc;
close all;

set(0,'defaultAxesFontSize',16);
addpath("GOCE_L1b_MatlabReaders_1.1/data/")
addpath("GOCE_L2b_MatlabReaders/data/")

% constants
Re = 6378.1370E3; % [m]
gpsEpoch = datetime(1980, 1, 6, 0, 0, 0, 'TimeZone', 'UTC');

date = "16_Nov_2012";

% load data
var = ["pos","vel","posCov","velCov", "time"];
for j = 1:5
file = var(j)+"_"+date+".mat";
load(file);
end
load(date + "_L2position.mat");

% extract the covariance
sigma2_pos = ones(4, length(TT_GPS_PVT_final));
sigma2_vel = ones(3, length(TT_GPS_PVT_final));

for j = 1:length(TT_GPS_PVT_final)
    % extract position uncertianty
    P = squeeze(COV_POS_FINAL(j, :, :));
    sigma2_pos(:, j) = diag(P);

     % extract velocity uncertianty
    P = squeeze(COV_VEL_FINAL(j, :, :));
    sigma2_vel(:, j) = diag(P);
end

% compute difference with level 2
[commonTimes, idx1, idx2] = intersect(TT_GPS_PVT_final, positions(:, 1));
diference = positions(idx2, 2:4) - POS_PVT_FINAL(idx1, :);    % [m]

% compute geographic coordinates
lla = ecef2lla(POS_PVT_FINAL);
lat = lla(:,1);
lon = lla(:,2);

% plot coverage
figure;
geoplot(lat, lon, '.b') 
geobasemap('landcover')
title('Sub-satellite Track')

% plot trajectory
figure()
subplot(1, 2, 1)
plot(TT_GPS_PVT_final, POS_PVT_FINAL./1E3, 'LineWidth', 2);
ylabel('[Km]')
xlabel('GPS sec')
title('Position ECEF')
legend('x', 'y', 'z')
grid on;

subplot(1, 2, 2)
plot(TT_GPS_PVT_final, VEL_PVT_FINAL./1E3, 'LineWidth', 2)
ylabel('[Km/s]')
xlabel('GPS sec')
title('Velocity ECEF')
legend('x', 'y', 'z')
grid on;
ylim([-10, 10])
sgtitle('GOCE orbit kinematic solution. Level 1 data')

% plot uncertainty
figure()
subplot(1, 2, 1)
plot(TT_GPS_PVT_final, sqrt(sigma2_pos(1:3, :)), 'LineWidth', 2)
ylabel('[m]')
xlabel('GPS sec')
title('Position ECEF')
legend('x', 'y', 'z')
grid on;

subplot(1, 2, 2)
plot(TT_GPS_PVT_final, sqrt(sigma2_vel), 'LineWidth', 2)
ylabel('[m/s]')
xlabel('GPS sec')
title('Velocity ECEF')
legend('x', 'y', 'z')
grid on;
sgtitle('GOCE orbit uncertainty. Level 1 data')

% plot positions difference between level 1 and level 2
time  = gpsEpoch + seconds(commonTimes);
figure()
plot(time, diference, 'LineWidth', 2)
ylabel('[m]')
xlabel('Time')
title('Position difference data level 1&2 ECEF')
legend('x', 'y', 'z')
grid on;

