clear;
clc;
close all;
format long g;
addpath('../functions/')

%%              EIG FRAME DEFINITION
% Description: try to define correctly the EIG frame given a gravity field
% with J2.

lambda = pi/20;
phi = pi/13;
theta = pi/15;
g = phi  + theta;
g2 = lambda + theta;

M3 = [cos(phi), sin(phi), 0;-sin(phi), cos(phi), 0;0, 0, 1];
M31= [cos(theta), sin(theta), 0;-sin(theta), cos(theta), 0;0, 0, 1]; 
M3t = [cos(g), sin(g), 0;-sin(g), cos(g), 0;0, 0, 1];

M2 = [cos(lambda), 0, -sin(lambda); 0, 1, 0;sin(lambda), 0, cos(lambda)];
M21 = [cos(theta), 0, -sin(theta); 0, 1, 0;sin(theta), 0, cos(theta)];
M2t = [cos(g2), 0, -sin(g2); 0, 1, 0;sin(g2), 0, cos(g2)];

L1 = M31 * M3 * M2;
L2 = M3t * M2;

Al = sym("Al","real");
th = sym("th","real");
Ap = sym("Ap","real");
G1 = sym("G1","real");
G12 = sym("G12","real");
G2 = sym("G2","real");
G3 = sym("G3","real");

Ar = pi/2;
M3 = [cos(Ap), sin(Ap), 0;-sin(Ap), cos(Ap), 0;0, 0, 1];
M3t = [cos(th), sin(th), 0;-sin(th), cos(th), 0;0, 0, 1];
M2 = [cos(Al), 0, -sin(Al); 0, 1, 0;sin(Al), 0, cos(Al)];
M2t = [cos(th), 0, -sin(th); 0, 1, 0;sin(th), 0, cos(th)];
M1 = [1, 0, 0;0, 0, 1;0, -1, 0];

L = M3 * M2 * M1;
T =[G1,0,0;0,G2,0;0,0,G3];
T = L * T * L';
T = simplify(T);

% Implement J2 & J3 effect

% position vector
r = 500;            % distance 2 planet [m]
l = deg2rad(75);     % longitude [rad]
p = deg2rad(40);    % latitude  [rad]

x = r * cos(p) * cos(l);
y = r * cos(p) * sin(l);
z = r * sin(p);
r_ECEF = [x;y;z];
v_ECEF = [0; 0.072111; 0.072111];

% planet params
GM = 5.2;           % point mass parameter [m^3/s^2]
Ref= 250;           % reference radius [m]
J2 = -0.4;          % J2 = -C20 coefficient
J3 = 0;           % J3 = -C30 coefficient
Cnm = [1, 0, 0, 0, 0, 0, 0;...
    0, 0, 0, 0, 0, 0, 0;...
    -J2, 0, 0, 0, 0, 0, 0;
    -J3, 0, 0, 0, 0, 0, 0];

% gravity gradient tensor components 3x3. Spherical coordinates
T_CL  = compute_GGT(r, p, l, GM, J2, J3, Ref);
[~, ~, T_ECEF] = potentialGradient_nm(Cnm, zeros(7,7), 3, ...
                                           r_ECEF, Ref, GM, 0);
[L, v] = eig(T_CL);
L = L';
theta  = atan2(L(3,2), L(2,2)); % L(2,2) or L(1,2)

g = theta;
M3t = [cos(g), sin(g), 0;-sin(g), cos(g), 0;0, 0, 1];
M2t = [cos(g), sin(g), 0; 0, 0, 1;-sin(g), cos(g), 0];

EIG_CL = M3t; % M3t or M2t

T_EIG = EIG_CL * T_CL * EIG_CL';

[CL_ECEF] = rotate2CL(p, l);
T_ECEF2 = CL_ECEF' * T_CL * CL_ECEF;

% compute RTN frame. Checking porpusses
[F, ~] = eig(T_ECEF);
r_ECEF = F(:, 3);
[ECEF_RTN] = RTN2ECI(r_ECEF, v_ECEF);
RTN_ECEF = ECEF_RTN';

% position vector
r2 = r + 10;               % distance 2 planet [m]
l2 = l + deg2rad(1E-1);    % longitude [rad]
p2 = p + deg2rad(1E-1);    % latitude  [rad]

T2_CL  = compute_GGT(r2, p2, l2, GM, J2, J3, Ref);

% rotate to ECEF
[CL_ECEF] = rotate2CL(p2, l2);
T2_ECEF = CL_ECEF' * T2_CL * CL_ECEF;

% rotate to T1 EIG frame
[CL_ECEF] = rotate2CL(p, l);
R = EIG_CL * CL_ECEF;
T2_EIG = R * T2_ECEF * R';

% compute sensitivity
[H_CL]  = compute_H(r, p, GM, Ref); 
H_ECEF = CL_ECEF' * H_CL * CL_ECEF;
H_EIG  = EIG_CL * H_CL * EIG_CL';
H_RTN  = RTN_ECEF * H_ECEF * ECEF_RTN; 

%%      FUNCTIONS
function [T]  = compute_GGT(r, p, l, GM, J2, J3, R) 
    % Description: return the grav. gradient tensor compoenents in
    % spherical coordinates including the J2 effect
    
    r3 = r^3;
    r5 = r^5;
    r6 = r^6;
    A = J2 * GM * R^2;
    B = J3 * GM * R^3;

    Trr = 2*GM/r3 - 6 * A / r5 * (3*sin(p)^2 - 1) + ...
        -10 * B / r6 * (5*sin(p)^3 - 3*sin(p));
    Trp = 6 * A /r5 * sin(2*p) + ...
        5/2 * B /r6 * (15*sin(p)*sin(2*p)/2 - 3*cos(p));
    Trl = 0;
    
    Tpp = -GM/r3 +3/2 * A /r5 * (7*sin(p)^2 - 3) + ...
        -1/2 * B /r6 * (3*sin(p) + 15*cos(p)*sin(2*p)/2 +30*sin(p)*cos(2*p)/2) + ...
        2 * B /r6 * (5*sin(p)^3 - 3*sin(p));
    Tpl = 0;

    Tll = -(Trr + Tpp);

    T = [Trr, Trp, Trl;Trp, Tpp, Tpl;Trl, Tpl, Tll];
end


function [T]  = compute_H(r, p, GM, R) 
    % Description: return the grav. gradient tensor compoenents in
    % spherical coordinates including the J2 effect
    
    r3 = r^3;
    r5 = r^5;
    r6 = r^6;
    A = 1 * GM * R^2; % J2
    B = 0 * GM * R^3; % J3
    C = 0;            % C00

    Trr = 2*GM/r3 * C - 6 * A / r5 * (3*sin(p)^2 - 1) + ...
        -10 * B / r6 * (5*sin(p)^3 - 3*sin(p));
    Trp = 6 * A /r5 * sin(2*p) + ...
        5/2 * B /r6 * (15*sin(p)*sin(2*p)/2 - 3*cos(p));
    Trl = 0;
    
    Tpp = -GM/r3 * C +3/2 * A /r5 * (7*sin(p)^2 - 3) + ...
        -1/2 * B /r6 * (3*sin(p) + 15*cos(p)*sin(2*p)/2 +30*sin(p)*cos(2*p)/2) + ...
        2 * B /r6 * (5*sin(p)^3 - 3*sin(p));
    Tpl = 0;

    Tll = -(Trr + Tpp);

    T = [Trr, Trp, Trl;Trp, Tpp, Tpl;Trl, Tpl, Tll];
end

function [CL_ECEF] = rotate2CL(p, l)
    Ar = pi/2;
    Ap = l;
    Al = p;
    M1 = [1, 0, 0;0, cos(Ar), sin(Ar);0, -sin(Ar), cos(Ar)];
    M3 = [cos(Al), sin(Al), 0;-sin(Al), cos(Al), 0;0, 0, 1];
    M2 = [cos(Ap), 0, -sin(Ap); 0, 1, 0;sin(Ap), 0, cos(Ap)];
    CL_ECEF = M3 * M2 * M1;
end
