clear; 
clc;
close all;

set(0,'defaultAxesFontSize',16);

addpath("data/")
addpath("functions/")
addpath("functions/solver")
addpath("functions/measurements")
addpath("functions/integrator")

% load SPICE kernels
cspice_furnsh('/Users/sergiocollibars/Documents/MATLAB/kernels/kernels.tm') 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       INVESTIGATE CORRELATIONS                          %
% Description: Invesitgate correlations along NRHO between states         %
% and when estimating extended states                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TBD