clear;
clc;
close all;

%% RESONANCES CALCULATOR
% Description: This code plots the semi-major values for N, M resonances.
% N: number of rotations on the primary
% M: number of orbits

% Input
N = [5, 6, 7, 8];
M = [1, 2];

N_elem = length(N)*length(M);

T = 4.297461 * 60 * 60;             % rotational period of the planet [s]
mu = 5.2;                           % grav. parameter [m^3/s^2]

% semi-major axis 
a_res = ones(N_elem, 1) * NaN;
label = {};

% loop
n = 1;
for i = 1:length(N)
    for j = 1:length(M)
        str = string(N(i)) + "/" + string(M(j));
        label{end+1} = str;
        a_res(n) = (N(i)/M(j) * T/(2*pi) * sqrt(mu))^(2/3);
        n= n + 1;
    end
end

% plot
figure();
plot(linspace(1, N_elem, N_elem), a_res./1000, 'Marker','diamond', ...
    'LineStyle', 'none', 'Color',"#D95319", 'MarkerFaceColor',"#D95319");
xlabel('Resonances');
ylabel('Semi- major axis, a [Km]');
title('Semi major axis for resonances at Bennu');
set(gca, 'XTick',1:N_elem, 'XTickLabel',label);
grid on;
