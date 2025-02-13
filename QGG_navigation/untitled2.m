clear;
clc;
close all;

cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm')

% Sample Data (Replace with actual data)

coords = load('postion.mat').posIterMC;  % Simulated random trajectory data
true_trajectory = load('true_postion.mat').state(:, 1:3)';  % Simulated true trajectory
time = load('time.mat').time;

% get scale parameters
[planetParams, ~, ~, ~, ~] = load_universe("CR3BP", ...
    [0, 1], 1/60);
cspice_kclear

T = length(coords(1, :, 1))/5;    % Number of time steps
N = length(coords(1, 1, :));    % Number of iterations

% dimensionalize vales:
coords = coords(:, :, :).*planetParams(2)/1E3;  % [km]
true_trajectory = true_trajectory.*planetParams(2)/1E3;  % [km]
time = time./planetParams(3)/3600;                      % [h]

% Video Setup
video_name = 'trajectory_animation.mp4';
v = VideoWriter(video_name, 'MPEG-4');
v.FrameRate = 60;
open(v);

% Plot Setup
figure;
hold on;
grid on;
xlabel('X [km]'); ylabel('Y [km]'); zlabel('Z [km]');
xlim([min(coords(1,:,:), [], 'all') max(coords(1,:,:), [], 'all')]);
ylim([min(coords(2,:,:), [], 'all') max(coords(2,:,:), [], 'all')]);
zlim([min(coords(3,:,:), [], 'all') max(coords(3,:,:), [], 'all')]);


% Plot Trajectories Over Time
for t = 1:T
    cla;  % Clear axes for each frame
    
    % Plot all iterations at current time step
    for n = 1:N
        plot3(coords(1,1:t,n), coords(2,1:t,n), coords(3,1:t,n), 'b', 'LineWidth', 1);
    end
    
    % Plot true trajectory
    plot3(true_trajectory(1,:), true_trajectory(2,:), true_trajectory(3, :), ...
          'r', 'LineWidth', 2);
    
    % Highlight current point
    scatter3(true_trajectory(1,t), true_trajectory(2,t), true_trajectory(3,t), ...
             70, 'k', 'filled');
    
    % Title with Time Step
    title(sprintf('Time : %.1f [hours]', time(t)));
    
    % Adjust View Angle (Optional)
    view([330 14]);  % Adjust as needed
    
    % Capture Frame
    frame = getframe(gcf);
    writeVideo(v, frame);
end

% Close Video
close(v);
disp('Video created successfully.');
