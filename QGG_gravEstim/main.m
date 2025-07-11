clear AttitudeData AttitudeObj estimData estimObj OrbitData OrbitObj AccObj;
clear file N name readData save SH_path SimData t_max t_min Time;
% % clc;        % WARNING: remember to comment before running MC simulation
close all;
format long g;

%%                                  MAIN
% ----------------------------------------------------------------------- %
%   Author: Sergio Coll Ibars
%
%   Date: 26/10/2022
%
%   Description: Main program. Calls all the defined modules and compute
%   the sensor acceleration.
% ------------------------------------------------------------------------%

% add paths
addpath('./modules/attitude_module/');
addpath('./modules/config_module/');
addpath('./modules/orbital_module/');
addpath('./modules/acc_module/');
addpath('./modules/estimation_module/');

addpath('./modules/attitude_module/functions');
addpath('./modules/config_module/functions');
addpath('./modules/orbital_module/functions');
addpath('./modules/acc_module/functions');
addpath('./modules/estimation_module/functions');
addpath('./modules/estimation_module/functions/solvers');
addpath('./modules/estimation_module/functions/filters');

addpath('data_files/');
addpath('config/');

% create module's objects
AttitudeObj = AttitudeClass;
OrbitObj = OrbitClass;
AccObj = AccClass;
estimObj = estimationClass;

% call config_module
file = "phaseA.txt";
SH_path = "HARMCOEFS_Bennu_OSIRIS_0.txt";
[SimData, OrbitData, AttitudeData, estimData] = config_module(file);


AttitudeObj.AttitudeData = AttitudeData;
OrbitObj.OrbitData = OrbitData;
OrbitObj.SH_path = SH_path;
estimObj.EstimData = estimData;

% extract simulation data
t_min = SimData(1);                         % min simulaiton time [s]
t_max = SimData(2);                         % max simulation time [s]

N = SimData(3);                             % # simulation points

AccObj.sigmaN_ii = SimData(4);              % acc STD noise diagonal [1/s^2]
AccObj.sigmaN_ij = SimData(5);              % acc STD noise out diagonal [1/s^2]
AccObj.noiseProfile = SimData(8);           % acc noise profile

% Define simulation time
Time = linspace(t_min, t_max, N);

if(SimData(6) == 0)
    readData = 0;
else
    readData = 1;
end

if(SimData(6) == 0)
    disp('run QGG measurements generator');

    % call attitude module
    save = true;
    AttitudeObj = attitude_module(AttitudeObj, Time, save);
    
    % call orbit module
    save = true;
    [OrbitObj, AttitudeObj] = orbit_module(OrbitObj, AttitudeObj, Time, ...
        save);

    % call acceleration module
    save = true;
    name = "accData";
    [AccObj] = acc_module(AttitudeObj, OrbitObj, AccObj, Time, save, ...
        name, readData);
end

% call estimation module if needed
save = true;
if(estimObj.EstimData(1) == 1)
    disp('run estimation');
    name = "estimData";
    estimObj = estimation_module(estimObj, OrbitObj, save, ...
        name, SimData(7));
end
disp('DONE!');

