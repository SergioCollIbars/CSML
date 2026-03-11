clear;
clc; close all;

% loading Kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm');

% Get initial and final time.
t0 = cspice_str2et({'2020-01-10 12:00:00 TDB'});
tf = cspice_str2et({'2020-01-10 04:00:00 TDB'});

TGT_BODY  = 'EARTH'; % ID 9:2 NRHO: -60000
time      = linspace(t0, tf, 1000); % this can be a vector too.
REF_FRAME = 'J2000';
SPEED_LIGHT_CORR = 'NONE';
OBSERVER  = 'SUN';

% pos & vel SPICE
[sc_SPICE, ~] = cspice_spkezr(TGT_BODY, time, REF_FRAME, ...
    SPEED_LIGHT_CORR, OBSERVER);
 
% clear kernels
cspice_kclear;