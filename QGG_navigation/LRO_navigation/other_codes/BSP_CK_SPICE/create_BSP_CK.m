%% CONVERT DATA FROM LRO SIM TO SPICE BSP/CK FILES
clear;
clc;
close all;

metaKernelPath = './kernels/kernels_files.tm';
cspice_furnsh(metaKernelPath);

% Add custom kernels for your simulated spacecraft
cspice_furnsh('lro_est.tsc');
cspice_furnsh('lro_est.tf');

%% Load MATLAB Sim data
load("./data/orbit_LRO_ideal_dataOut.mat");

estim_state = mC_struct.sim1.Xf;
true_state  = mC_struct.sim1.state_true;
t_et        = time;

BN_mat = compute_orientation_SC(time, estim_state(1:6, :)', ...
    mC_struct.sim1.I_ALIG, NB_MOON_mat);

%% Create BSPs
SC_EST_ID  = -999000;
SC_TRUE_ID = -999001;

writeTrajectoryToBSP(t_et, estim_state(1:3, :)./1E3, ...
    estim_state(4:6, :)./1E3, 'estim_traj.bsp', SC_EST_ID);

writeTrajectoryToBSP(t_et, true_state(1:3, :)./1E3, ...
    true_state(4:6, :)./1E3, 'true_traj.bsp', SC_TRUE_ID);

%% Create CK for estimated spacecraft attitude
write_dcm_ck(t_et, BN_mat, 'estim_attitude.bc', SC_EST_ID, -999010, 'J2000');

cspice_kclear();