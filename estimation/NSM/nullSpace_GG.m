clear;
clc;
close all;
addpath('../functions/');


% generate tensor meas. 
n_max = 0;
normalized = 0;
Re = 250;
GM = 5.2;

% distance vector. ECEF (truth)
r = [1000;100;50];
Ar = 1E-3.*[1;1;1];     % [m]
At = 0.*[1;1;1];   % [rad]
DGM = 0.0005;

Cmat = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0.1 0.2 0 0 0 0;...
        0.2 0.3 0.04 0 0 0 0;...
        0.02 0.0001 -0.0004 -7-05 4.51-05 0 0;...
        -0.005 0.00044 9-05 1e-05 1e-05 1e-06 0;...
        -0.009 0.001 9e-05 2e-05 -2e-06 9e-08 2e-07];

Cmatap = [1 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        -0.4 0.1 0.2 0 0 0 0;...
        0.2 0.3 0.04 0 0 0 0;...
        0.02 0.0001 -0.0004 -7-05 4.51-05 0 0;...
        -0.005 0.00044 9-05 1e-05 1e-05 1e-06 0;...
        -0.009 0.001 9e-05 2e-05 -2e-06 9e-08 2e-07];
Ac = zeros(45, 1);
Ac(1) = 0.001;

Smat = [0 0 0 0 0 0 0;...
        0 0 0 0 0 0 0;...
        0 1E-6 -3e-05 0 0 0 0;...
        0 0.0001 7e-05 -0.001 0 0 0;
        0 0.007 0.0001 0.0001 0.003 0 0;...
        0 -2e-05 -8e-05 5e-06 -0.0005 0.0005 0;...
        0 -1e-06 -1e-05 -8e-06 -2e-06 -0.0006 -4e-05];

% measurements
sigma = 0;
[~, ~, ddU1] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                r, Re, GM, ...
                                                normalized);
ddU1 = ddU1 + normrnd(0, sigma, [3, 3]);
% rotate measurements
[R] = rotationMatrix(At(1), At(2), At(3), [1, 2, 3]);
ddU1 = R' * ddU1 * R;

[~, ~, ddU2] = potentialGradient_nm(Cmatap, Smat, n_max, ...
                                                r+Ar, Re, GM+DGM, ...
                                                normalized);
% H partials w.r.t position
[Hpos] = compute_posPartials(n_max, normalized, Cmatap, Smat, Re, GM+DGM, r+Ar);
posVec = r+Ar;
j =ones(6, 1) * Ar(1)^2;
[Hpos2] = compute_posPartials_2ndOrder(GM, posVec(1), posVec(2), posVec(3));

% H partials w.r.t SH coeff. Using 4 meas only
[~, Hc] = potentialGradient_Cnm(n_max, r+Ar, Re, 1, eye(3,3), 0); 

% H partials w.r.t attitude
[Hrot] = compute_rotPartials(n_max, normalized, Cmatap, Smat, Re, GM+DGM, r+Ar);

% meas difference
Y = ddU1 - ddU2;
dY1 = [Y(1,1);Y(1,2);Y(1,3);Y(2,2);Y(2,3)];
% % dY = [Y(1,1);Y(2,2);Y(3,3)];

% compute null Space
c = [Hpos(1, :);Hpos(2, :);Hpos(3, :);Hpos(5, :);Hpos(6, :)];
d = [Hrot(1, :);Hrot(2, :);Hrot(3, :);Hrot(5, :);Hrot(6, :)];
b = [Hc(1, 1:end); Hc(2, 1:end); Hc(3, 1:end); Hc(5, 1:end); Hc(6, 1:end)];
C = null(c');          % H' * C = 0 // C' * H = 0
B = null(b');
D = null(d');  

% stack null space
Nspace = 30;     % number of time steps
dr  = linspace(10, 100, Nspace);
Bs = ones(length(dr)*5, length(Hc) - 1);
Cs = ones(length(dr)*5, 3);
Ds = ones(length(dr)*5, 3);
dYs = ones(length(dr)*5, 1);
for j = 1:length(dr)
    dR = [dr(j);dr(j);dr(j)];

    % min/max index
    max = 5 * j;
    min = max - 4;

    % compute meas. increment
    [~, ~, ddU1] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                r+dR, Re, GM, ...
                                                normalized);
    ddU1 = ddU1 + normrnd(0, sigma, [3, 3]);
    [R] = rotationMatrix(At(1), At(2), At(3), [1, 2, 3]);
    ddU1 = R' * ddU1 * R;

    [~, ~, ddU2] = potentialGradient_nm(Cmatap, Smat, n_max, ...
                                                    r+Ar+dR, Re, GM + DGM, ...
                                                    normalized);
    Y = ddU1 - ddU2;
    dYs(min:max) = [Y(1,1);Y(1,2);Y(1,3);Y(2,2);Y(2,3)];

    % compute & stack position partials
    [Hpos] = compute_posPartials(n_max, normalized, Cmatap, Smat, Re,...
        GM + DGM, r+Ar+dR);

    c = [Hpos(1, :);Hpos(2, :);Hpos(3, :);Hpos(5, :);Hpos(6, :)];

    % H partials w.r.t attitude
    [Hrot] = compute_rotPartials(n_max, normalized, Cmatap, Smat, Re, GM, r+Ar+dR);
    d = [Hrot(1, :);Hrot(2, :);Hrot(3, :);Hrot(5, :);Hrot(6, :)];
    Ds(min:max, :) = d;
    Cs(min:max, :) = c;

    % compute & stack gravity field partials
    [~, H] = potentialGradient_Cnm(n_max, r+Ar+dR, Re, GM + DGM, eye(3,3), 0); 

    Bs(min:max, :) = [H(1, 2:end); H(2, 2:end); H(3, 2:end); ...
        H(5, 2:end); H(6, 2:end)];
end

% compute stack null spaces
C = null(Cs');  % position & attitude null space
B = null(Bs');  % gravity field null space
D = null(Ds');  % attitude null space
%% FUNCTIONS
function [Hpos] = compute_posPartials(n_max, normalized, Cmat, Smat, Re, GM, r)
    % output value
    Hpos = ones(9, 3) * NaN;
    
    eps = 1E-5;
    for j = 1:3
        Ar = zeros(3, 1);
        Ar(j) = eps;

        rpos = r + Ar./2;
        rneg = r - Ar./2; 

        [~, ~, ddUpos] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                rpos, Re, GM, ...
                                                normalized);
        [~, ~, ddUneg] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                rneg, Re, GM, ...
                                                normalized);
        H = (ddUpos - ddUneg)./(eps);
        Hpos(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
            H(3, 1); H(3,2); H(3,3)];
    end
end

function [Hrot] = compute_rotPartials(n_max, normalized, Cmat, Smat, Re, GM, r)
    % output value
    Hrot = ones(9, 3) * NaN;
    
    eps = 1E-5;
    for j = 1:3
        At = zeros(3, 1);
        At(j) = eps;

        Atpos = At./2;
        Atneg = - At./2; 

        [Rpos] = rotationMatrix(Atpos(1), Atpos(2), Atpos(3), [1, 2, 3]);
        [Rneg] = rotationMatrix(Atneg(1), Atneg(2), Atneg(3), [1, 2, 3]);

        [~, ~, ddUpos] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                r, Re, GM, ...
                                                normalized);

        ddUpos = Rpos' * ddUpos * Rpos;

        [~, ~, ddUneg] = potentialGradient_nm(Cmat, Smat, n_max, ...
                                                r, Re, GM, ...
                                                normalized);
        ddUneg = Rneg' * ddUneg * Rneg;

        H = (ddUpos - ddUneg)./vecnorm(Atpos-Atneg);
        Hrot(:, j) = [H(1,1);H(1,2);H(1,3);H(2,1);H(2,2);H(2,3);...
            H(3, 1); H(3,2); H(3,3)];
    end
end