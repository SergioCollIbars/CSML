clear;
clc;
close all;

%% PLANET EFFECT IN THE CR3BP

% Unit conversion
AU =  1.495978707E11;       % [m]
Etv = 1E-9;                 % [1/s^2]
D = 384399e3;               % [m]
GM1 = 398600.44150E9;       % [m^3/s^2]
GM2 = 4902.80011526323E9;   % [m^3/s^2]
mu = GM2 / (GM1 + GM2);     % [-]
Earth_BC = mu*D;        % [m]

% planets distance to Earth
planets = ["Mercury", "Venus" "Mars", "Sun", "Jupiter", " Saturn", "Uranus", "Neptune", "Pluto"];
PD = [0.552; 0.266; 0.372; 0.983; 3.957; 8.050; 17.292; 28.817; 28.699 ];    % [AU]

% planets GM [m^3/s^2]
GMP = [22031.86855;324858.592;42828.37582;132712E6;126712764.1;37940584.84;5794556.4;6836527.101;975.5].*1E9;

FNP = zeros(1, length(PD));
for j = 1:length(PD)
    x = (PD(j)*AU + Earth_BC);% [m]
    y = 0; 
    z = 0;                     % [m]
    GM = GMP(j);               % [m^3/s^2]

    [ddU] = compute_QGG(GM, x, y, z);   % [1/s^2]
    ddU = ddU./Etv;                     % [E]

    [FNP(j)] = compute_FrobeniusNorm(ddU);
end

figure()
bar(FNP)
set(gca,'xticklabel',planets)
set(gca,'XScale','log')
ylabel('Eotvos [E]');
title('Frobenius norm at CR3BP barycenter')

%% FUNCTIONS
function [ddU] = compute_QGG(GM, x, y, z)
    rvec = [x;y;z];
    r3 = vecnorm(rvec)^3;
    r5 = vecnorm(rvec)^5;

    U11 = -3*x*x/r5 + 1/r3;
    U12 = -3*x*y/r5;
    U13 = -3*x*z/r5;
    U22 = -3*y*y/r5 + 1/r3;
    U23 = -3*y*z/r5;
    U33 = -3*z*z/r5 + 1/r3;

    ddU = GM * [U11, U12, U13;...
                U12, U22, U23;...
                U13, U23, U33];
end

function [FN] = compute_FrobeniusNorm(A)
    FN = (trace(A'*A))^(0.5);
end
