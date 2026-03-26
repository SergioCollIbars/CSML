import numpy as np
from functions.read_GG_obs import read_GG_obs
from functions.GG_rotation_partials import GG_rotation_partials
"""
TEST NULL SPACE COMPUTED THROUGH SVD (in PYTHON)

Author: Sergio Coll-Ibars
Date:   05/03/2026

Description: Use the obervation files to compute the 
        orientation sensitivity matrix for the GG.
        The compute the null space through SVD.
        If the MATLAB file exist already, do comparison.
        Otherwise, just save the partials in a txt file.
"""

INPUT_OBS = "/Users/sergiocollibars/Documents/GG_observations/120by120/" + \
            "goce_eggreg_2008-08-01_RL5061_RL05surpv6.980.001.ggr"

MATLAB_NS = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/modify_regres_codes/MATLAB_version/"+\
            "NS_MATLAB_2008_08_01.txt"
matlab_data = np.loadtxt(MATLAB_NS)

G = read_GG_obs(INPUT_OBS)

Nt  = len(G)
Nobs = len(G.T)

NS_matlab = matlab_data.reshape((Nt, 6, 3))

for k in range(Nt):
    g_obs = G[k, 1:Nobs]
    Hrot = GG_rotation_partials(g_obs)

    # SVD for the H_rot
    U, S, Vh = np.linalg.svd(Hrot.T)
    V  = Vh.T # NOTE: Pyhton returns the Hermitian. We need to transpose

    P  = V[:, 3:6]
    NS = NS_matlab[k] 

    diff = np.max((P - NS) / NS)
    print(f"Max rel. error MATLAB vs PYTHON null space. Epoch {G[k, 0]} is {diff}")

print('Execution finished!')