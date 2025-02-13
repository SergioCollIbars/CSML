clear;
clc;
format long g;
%%                  EULER TRIANGLE METHOD

% theta angles
th1 = deg2rad(2.017); % [rad] Omega
th2 = deg2rad(6.034); % [rad] i
th3 = deg2rad(66.304 + 64.541); % [rad] omega + f

% Phi angles
ph1 = deg2rad(-90); % [rad] Omega
ph2 = deg2rad(90); % [rad] i
ph3 = deg2rad(90); % [rad] omega

% theta angles
th1 = deg2rad(0); % [rad] Omega
th2 = deg2rad(10.8); % [rad] i
th3 = deg2rad(0); % [rad] omega + f

% Phi angles
ph1 = deg2rad(0); % [rad] Omega
ph2 = deg2rad(90); % [rad] i
ph3 = deg2rad(0); % [rad] omega

% 3-1-3 addition (th1, th2, th3) + (ph1, ph2, ph3) = (ps1, ps2, ps3)
ps2 = acos(cos(th2)*cos(ph2) - sin(th2)*sin(ph2)*cos(th3+ph1));

a = sin(th2)*sin(ph2)*sin(th3+ph1);
b = cos(ph2) - cos(th2)*cos(ps2);
ps1 = th1 + atan(a/b);

a = sin(th2)*sin(ph2)*sin(th3+ph1);
b = cos(th2) - cos(ph2)*cos(ps2);
ps3 = ph3 + atan(a/b);

% disp values
disp("Psi 1 (Omega): " + string(rad2deg(ps1)));
disp("Psi 2 (i): " + string(rad2deg(ps2)));
disp("Psi 3 (omega): " + string(rad2deg(ps3)));

SAR_N = [-0.680087366936192, 0.728805664194894 ,0.0795202940841892;...
     -0.733121739641444, -0.676613679072915, -0.0687491393147374;...
      0.00369975660015884, -0.105053477471489,  0.994459691828808];
% 
% SAR_N = [0.800139124577202, 0.596799112578933,0.0600682990089447;...
%         -0.599803045275638,0.795483869161497,0.086265408953561;...
%        0.00369975660015885,-0.105053477471489,0.994459691828808];


B_SAR = [0 ,0 ,1;...
     -0, 1, 0;...
      -1, 0, 0];

BN = B_SAR * SAR_N;

Omega = atan2(BN(3,1), -BN(3,2));
omega = atan2(BN(1,3), BN(2,3));
i = acos(BN(3,3));

disp("Omega: " + string(rad2deg(Omega)));
disp("i: " + string(rad2deg(i)));
disp("omega: " + string(rad2deg(omega)));