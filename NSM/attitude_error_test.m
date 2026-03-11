clear; clc; close all;

%% Settings
fs = 1/5;         % Hz
T  = 2000000;      % seconds

plotCommand = 1;    % 1: yes / 0: no
[x] = create_noise_from_PSD(fs, T, plotCommand);