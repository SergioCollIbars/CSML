clear;
clc;
close all;

%%                      PLOT COVARIANCE 
% Description: Plot the covariance againts the MC simulation resutls

% PATH
addpath("functions/")
addpath("//Users/sergiocollibars/Desktop/CSML/codes/QGG/data_files/")
addpath("//Users/sergiocollibars/Desktop/GOCE_analysis/phaseB/obtained/LS/")
addpath("//Users/sergiocollibars/Desktop/GOCE_analysis/phaseB/requiered/LS/")

% MODE
mode = 2;                       % obtained values = 1 / required values = 2

% LOAD DATA
if(mode == 1)
    err = load("obtained_error_MC.mat").CoefErr;
    errRMS = load("obtained_RMS_error_MC.mat").CoefErr_RMS;
elseif(mode == 2)
    err = load("required_error_MC.mat").CoefErr;
    errRMS = load("required_RMS_error_MC.mat").CoefErr_RMS;
end
Ncs = length(err(1, :));

T4 = readtable("estimData_4.txt").P;
T6 = readtable("estimData_6.txt");
Px = reshape(T4, [58, 58]);

% compute statistics
Px = Px(1:Ncs, 1:Ncs);
sigma = sqrt(diag(Px));
[sigmaRMS] = computeRMS_coefErr(6, 26, 20, sigma);

bound_up = 1*sigmaRMS + median(errRMS); 
bound_down = -1*sigmaRMS + median(errRMS); 

figure()
plot(2:6, errRMS(:, 2:end), 'Marker','o', ...
    'LineStyle', 'none', 'Color','r', 'MarkerFaceColor','r')
hold all;
plot(2:6 ,bound_up(2:end), 2:6 ,bound_down(2:end),...
    'LineWidth', 1.5, 'Color', 'g')
grid on;
title('RMS error compared with truth. MC simualtion = 200')
xlabel('Harmonics degree [n]')
ylabel('coefficient error')
xlim([1.75, 6.25])

% Phase B
%   obtained : 6.32E-11 / 2.52E-9
%   required : 5.06E-12 / 5.06E-12

