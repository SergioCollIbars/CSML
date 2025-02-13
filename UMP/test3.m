clear;
clc;
close all;

% load data files
addpath('//Users/sergiocollibars/Desktop/data');


% GM TEST
truth = load('truthVal.mat').RSS_sigma';

data1 = load('deltaGM_val.mat').RSS_sigma';
data2 = load('deltaGM_val_2.mat').RSS_sigma';
deltaGM = [data1, data2];

deltaGM_error = zeros(length(truth(:, 1)), length(deltaGM(1, :)));
for col = 1:length(deltaGM(1, :))
    deltaGM_error(:, col) = abs((truth - deltaGM(:, col))./(truth));
end

openfig('deltaGM_error.fig');
hold on;
x = [0.2, 0.5];
plot(x, deltaGM_error(1, :), 'o', 'Color', 'b');
legend('Anlytical', 'Numerical')


% Eccentrycity test
truth = load('truth.mat').RSS_sigma';
data1 = load('e_val.mat').RSS_sigma';
data2 = load('e_val_2.mat').RSS_sigma';
data3 = load('e_val_3.mat').RSS_sigma';
data4 = load('e_val_4.mat').RSS_sigma';
data5 = load('e_val_5.mat').RSS_sigma';
e_val = [data1, data2, data3, data4, data5];

e_error = zeros(length(truth(:, 1)), length(e_val(1, :)));
for col = 1:length(e_val(1, :))
    e_error(:, col) = abs((truth - e_val(:, col))./(truth));
end

openfig('e_error.fig');
hold on;
x = [0.05, 0.1, 0.02, 0.03, 0.15];
semilogy(x, e_val(:, :), 'o', 'MarkerFaceColor','auto');
legend('Anlytical', 'Numerical')
xlim([0, 0.25])

% Inclination test
truth = load('truthVal.mat').RSS_sigma';

data1 = load('i_val_5.mat').sigma_RSS;
data2 = load('i_val_10.mat').sigma_RSS;
data3 = load('i_val_20.mat').sigma_RSS;
data4 = load('i_val_50.mat').sigma_RSS;
i_val = [data1, data2, data3, data4];

i_error = zeros(length(truth(:, 1)), length(i_val(1, :)));
for col = 1:length(i_val(1, :))
    i_error(:, col) = abs((truth - i_val(:, col))./(truth));
end

openfig('i_error.fig');
hold on;
x = [5, 10, 20, 50];
plot(x, i_val(:, :), 'o');
legend('Anlytical', 'Numerical')
