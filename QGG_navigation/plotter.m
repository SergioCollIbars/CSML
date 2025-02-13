clear;
clc;
close all;

set(0,'defaultAxesFontSize',16);

% Load data
timeVector = load('time2.mat').time; % [Time]
positionMatrix = load('traj.mat').pos; % [Km]
m = load('measurement.mat').T; % [Eotvos]
numPoints = length(timeVector);

% format data
meas = ones(1, length(timeVector)) * NaN;
for j = 1:length(timeVector)
    v =  [m(1, j), m(2, j), m(3, j);...
        m(2, j), m(4, j), m(5, j);...
        m(3, j), m(5, j), m(6, j)];
    meas(j) = norm(v, 'fro');
end
time = timeVector;

% Compute vector norms
vecNorms = vecnorm(positionMatrix);

% Initialize figure
figure;
set(gcf, 'Color', 'w'); % Set background to white

% Video settings
movFilename = 'evolution';
video = VideoWriter(movFilename, 'MPEG-4'); % Create a video object
video.FrameRate = 936; % Set frame rate
open(video); % Open the video file for writing

% Loop through each point and update plots
for i = 1:numPoints
    % First subplot: Evolution of time vector
    subplot(2, 1, 1);
    plot(time(1:i), meas(1:i), 'g-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    title('Frobenius norm for gradiometer signal');
    xlabel('Time');
    ylabel('[Eotvos]');
    grid on;
    xlim([min(time), max(time)]);
    ylim([min(meas), max(meas(1:i))]);
    
    % Second subplot: Evolution of vector norm
    subplot(2, 1, 2);
    plot(time(1:i), vecNorms(1:i), 'b-o', 'LineWidth', 1.5, 'MarkerSize', 4);
    title('Relative postion w.r.t the Moon');
    xlabel('Time');
    ylabel('[km]');
    grid on;
    xlim([min(time), max(time)]);
    ylim([0, max(vecNorms)]);
    
    % Capture the frame and write to video
    frame = getframe(gcf);
    writeVideo(video, frame);
end

% Close the video file
close(video);

disp(['Video saved as ', movFilename]);
