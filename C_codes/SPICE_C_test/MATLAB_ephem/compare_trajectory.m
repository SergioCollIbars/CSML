clear;
clc;

% read CSPICE integrated data
data = readmatrix('MATLAB_ephem/trajectory_output.txt');
timeC = data(:, 1);
posC  = data(:, 2:4); % [m]
velC  = data(:, 5:7); % [m/s]

% Load matlab values
load("MATLAB_ephem/state.mat"); load("MATLAB_ephem/time.mat");

% compute errors
timeErr = time - timeC;

posErr  = state(:, 1:3) - posC; 
velErr  = state(:, 4:6) - velC; 

figure(); tt = time - time(1);
plot(tt, timeErr); title('Time error');

figure();
plot(tt, posErr); title('Position error');

figure();
plot(tt, velErr); title('Velocity error');