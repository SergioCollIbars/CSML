clear;
clc;
close all;

addpath("/Users/sergiocollibars/Desktop/CSML/MSODP_functions/att_errors_res_NSM/functions");

%% CREATE NULL SPACE MAPPING
% Description: using the GG observation for a particular day, compute and
% save the NS using SVD in a txt file. 
% This is use as reference for comparison with other codes.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Output file
output_file = "NS_MATLAB_2008_08_01.txt";

% Observation file
GG_obs_file = "/Users/sergiocollibars/Documents/GG_observations/120by120/" + ...
    "goce_eggreg_2008-08-01_RL5061_RL05surpv6.980.001.ggr";

mask = [1, 1, 1, 1, 1, 1];
% Read observation data
G = read_GG_obs(GG_obs_file);

Nt     = length(G(:, 1));
NS_val = zeros(Nt * 6, 3); 
for k = 1:Nt
    g = G(k, 2:end);
    g_T = [g(1);g(2);g(3);g(2);g(4);g(5);g(3);g(5);g(6)];

    [Hrot] = compute_rotPartials_analy(g_T, eye(3));

    h = [Hrot(1:3, :);Hrot(5:6, :);Hrot(9, :)];
    h = h(logical(mask), :);

    [s,v,d] = svd(h');

    V = d(:, 4:end);

    P = eye(sum(mask)) - h * ((h'*h) \ h');

    maxInd = 6 * k; minInd = maxInd - 5;
    NS_val(minInd:maxInd, :) = V;
end

% Save data
writematrix(NS_val, output_file, 'Delimiter', 'space');