clear;
clc;
close all;
format long g;


%%      EVALUATE NOISE PROJECTION
% Description: Undertand and evaluate noise projection in Null space
% method.
% Author: Sergio Coll
% Date: 10/10/24

% Asteroid parameters.
path = "HARMCOEFS_BENNU_OSIRIS_1.txt";
name = "BENNU";
[Cnm, Snm, Re] = readCoeff(path);
GM = 5.2;
n_max  = 6;
normalized = 1;
W = 4.06130329511851E-4;  % Rotation ang. vel   [rad/s]
W0 = 0;                   % Initial asteroid longitude
RA = deg2rad(86.6388);    % Right Ascension     [rad]
DEC = deg2rad(-65.1086);  % Declination         [rad]


poleParams = [W, W0, RA, DEC];
asterParams = [GM, Re, n_max, normalized];

% SH harmonics
[Nc, Ns, Ncs] = count_num_coeff(n_max); 

% Initial conditions
r      = 0.3E3;
phi    = pi/2;
lambda = 0;
theta  = pi/2 - phi;% Orbit colatitude [m]
R = [sin(theta)*cos(lambda), cos(theta)*cos(lambda), -sin(lambda);...
    sin(theta)*sin(lambda), cos(theta)*sin(lambda), cos(lambda);...
    cos(theta), -sin(theta), 0];
r0 = R * [r;0;0];           % [ACI]
v0 = R * [0;0;sqrt(GM/r)];  % [ACI]

% position error
Ar = 0.5*[1;1;1];            % [ACI]

% time vector
n = sqrt(GM / r^3);    % Mean motion         [rad/s]
T = (2 * pi / n);
rev = 2;
f = 1/60;
t = linspace(0, rev*T, rev*T * f);
Nt = length(t);

% integrate trajectory
options = odeset('RelTol',1e-13,'AbsTol',1e-13);
[~, state_t] = ode113(@(t, x) EoM(t, x, Cnm, Snm, n_max, GM, Re, normalized, ...
    W0, W, RA, DEC), t, [r0;v0], options);
rn = state_t(:, 1:3)';
vn = state_t(:, 4:6)';

tt = 'Orbit radius along trajectory. T = ' + string(T./3600) + ' h';
plot_orbit(state_t, 'BENNU', t./T ,Re, tt)


% create random points
% noise values
Nr = 1E4;
noise = zeros(5, Nr, Nt);
sigma   = [1E-12 1E-11, 1E-12, 1E-11, 1E-13];

sigma1  = 0.01 * 1E-9 * sqrt(f); % Vxx, Vyy
sigma2  = 0.6  * 1E-9 * sqrt(f); % Vyz, Vyx
sigma3  = 0.02 * 1E-9 * sqrt(f); % Vxz, Vzz

R = diag([sigma1, sigma2, sigma3, sigma1, sigma2].^2);

means    = zeros(1, 5);
std_devs = [sigma1, sigma2, sigma3, sigma1, sigma2]; 
num_realizations = Nr; % Number of realizations

x_corr  = zeros(Nr, Nt); y_corr  = x_corr;
x_dcorr = zeros(Nr, Nt); y_dcorr = x_dcorr;
for k = 1:length(t)
    noise(:, :, k) =  normrnd(repmat(means', 1, num_realizations), ...
    repmat(std_devs', 1, num_realizations));

    % position partials & null space
    [Hp] = compute_posPartials(n_max, normalized, Cnm, Snm, Re, GM, rn(:, k), ...
         eye(3,3));
    C = null([Hp(1, :);Hp(2,:);Hp(3,:);Hp(5, :);Hp(6, :)]');
    
    % noise projection
    r = C' * R * C;
    np = C' * noise(:, :, k);
    x_corr(:, k) = np(1, :);
    y_corr(:, k) = np(2, :);

    % noise decorrelated
    [v, ~] = eig(r);
    rd = v'*r*v;
    nd = v'*np;
    x_dcorr(:, k) = nd(1, :);
    y_dcorr(:, k) = nd(2, :);
end

% do-video
createVideo(x_corr, y_corr, x_dcorr, y_dcorr, t, vecnorm(rn), Nt, 3*sigma2);

% scatter plot
figure;
subplot(1, 2, 1)
scatter(np(1,:), np(2,:), 'filled');
title('Correlated');
xlabel('X-axis');
ylabel('Y-axis');
grid on;
axis equal;  % Equal scaling on both axes

subplot(1, 2, 2)
scatter(nd(1,:), nd(2,:), 'filled');
title('De-correlated')
xlabel('X-axis');
ylabel('Y-axis');
grid on;
axis equal;  % Equal scaling on both axes

sgtitle('2D Gaussian Cloud of Random Points');


%% FUNCTION

function [] = plot_orbit(state_t, name, time, Re, tt)
    plot_trajectory(state_t, name);
    Nt = length(time);
    figure()
    plot(time, vecnorm(state_t(:, 1:3)'), 'LineWidth', 2)
    hold on;
    plot(time, ones(1, Nt)*Re, 'LineWidth', 2, 'Color', 'r', 'LineStyle','--')
    xlabel('Orb. Period number, T')
    ylabel('[m]')
    title(tt)
    legend('orbit radius', 'brillouin sphere')
end

function [] = createVideo(x1, y1, x2, y2, time, pos, N,axlim)
    % Parameters for video
    filename = 'plot_evolution_video.mp4';  % Video file name
    fps = 10;  % Frames per second
    
    % Create video writer object
    video = VideoWriter(filename, 'MPEG-4');  
    video.FrameRate = fps;  % Set frame rate
    open(video);  % Open video file for writing
    
    % Create a figure for the plot
    figure;
    hold on;

    % Loop to create and capture frames
    for t = 1:N
        subplot(2, 2, 1)
        scatter(x1(:,t), y1(:,t), 'filled');
        title('Correlated');
        xlabel('V_1 - axis');
        ylabel('V_2 - axis');
        grid on;
        axis equal;  % Equal scaling on both axes
        axis([-axlim axlim -axlim axlim]);
        
        subplot(2, 2, 2)
        scatter(x2(:,t), y2(:,t), 'filled');
        title('De-correlated')
        xlabel('V_1 - axis');
        ylabel('V_2 - axis');
        grid on;
        axis equal;  % Equal scaling on both axes
        axis([-axlim axlim -axlim axlim]);

        sgtitle(sprintf('Frame: %d', t));  % Add a frame number to the title
        grid on;

        subplot(2, 2, [3 4])
        plot(time(1:t), pos(1:t), 'b-', 'LineWidth', 2);
        axis([time(1) time(end) 250 500]);
        hold on;  % Keep the plot for each frame
        xlabel('TIME [sec]')
        ylabel('[m]')
        title('Orbit radius')

        % Capture the current frame
        frame = getframe(gcf);  
        writeVideo(video, frame);  % Write the frame to the video
        
        pause(0.05);  % Small pause to simulate real-time plotting (optional)
    end
    
    % Close the video file
    close(video);
    
    disp('Video saved successfully!');
end
