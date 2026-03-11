import numpy as np
from functions.read_GG_obs import read_GG_obs

"""
TEST READ OBSERVATION FUNCTION

Author: Sergio Coll-Ibars
Date:   05/03/2026

Description: Compate the read observation function
    to the standard Python loadtxt() built in library
    Display in console the max. relative error between 
    both read observations per epoch 
"""

INPUT_OBS = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/GG_observations/" + \
            "goce_eggreg_2008-08-01_RL5061_RL05surpv6.980.001.ggr"


obs_data = np.loadtxt(INPUT_OBS)

ref_data = obs_data[:, 20:26]

G = read_GG_obs(INPUT_OBS)

diff = G[:, 1:len(G.T)] - ref_data
rel_diff = diff / ref_data

print('Relative error in the read observation function')
Nt = len(G)
for i in range(Nt):
    print(f"Max error at epoch {i} : {np.max(rel_diff[i, :])}")


