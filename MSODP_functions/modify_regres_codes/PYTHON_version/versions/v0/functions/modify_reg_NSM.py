import numpy as np
from functions.read_reg_file_partials import read_reg
from functions.read_GG_obs import read_GG_obs
from functions.GG_rotation_partials import GG_rotation_partials
from functions.update_reg_partials_allRows import reg_update_partials
import os

def modify_reg_NSM(GG_obs, reg_files, outDir):
    """
    from the regress files for the XX, XY, XZ, YY, YZ, ZZ components and 
    the GG observations, construct the Null Space and modify the partials
    in the regres fields.

    Save the new regress files in the output folder location
    name: 

    <original_XX_name.reg> --> <original_NSM_i_name.reg>

    Returns
    ----------
    N/A
    """

    # extract and load the XX, XY, XZ, YY, YZ, ZZ partials
    H_XX = read_reg(reg_files[0])
    H_XY = read_reg(reg_files[1])
    H_XZ = read_reg(reg_files[2])
    H_YY = read_reg(reg_files[3])
    H_YZ = read_reg(reg_files[4])
    H_ZZ = read_reg(reg_files[5])

    # extract and load GG observations
    G    = read_GG_obs(GG_obs)

    Nt        = len(H_XX)
    Nobs      = len(G.T)
    Npart     = len(H_XX.T)
    
    # projected partials
    H_1_proj = np.zeros((Nt, Npart))
    H_2_proj = np.zeros((Nt, Npart))
    H_3_proj = np.zeros((Nt, Npart))

    # project partials
    print("Computing NSM projection to regres files")
    for k in range(Nt):
        # partials
        h = np.array([H_XX[k, 0:Npart-2],H_XY[k, 0:Npart-2],H_XZ[k, 0:Npart-2],\
                    H_YY[k, 0:Npart-2],H_YZ[k, 0:Npart-2],H_ZZ[k, 0:Npart-2]])
        # residuals
        r = np.array([H_XX[k, Npart-2],H_XY[k, Npart-2],H_XZ[k, Npart-2],\
                    H_YY[k, Npart-2],H_YZ[k, Npart-2],H_ZZ[k, Npart-2] ]).reshape(6,1)
        # noise
        n = np.diag([H_XX[k, Npart-1], H_XY[k, Npart-1], H_XZ[k, Npart-1], \
                     H_YY[k, Npart-1], H_YZ[k, Npart-1], H_ZZ[k, Npart-1]])
        # observations
        g = G[k, 1:Nobs]

        # compute rotation partials
        Hrot = GG_rotation_partials(g)

        # SVD for the H_rot
        U, S, Vh = np.linalg.svd(Hrot.T)
        V = Vh.T # NOTE: Pyhton returns the Hermitian. (need to transpose for left NS)

        # contstrunc projector (Orthogonal version)
        P = V[:, 3:6]
        
        # project partials, residuals & sigmas
        h_p  = P.T @  h
        r_p  = P.T @  r
        n_ax = n @ P
        n_p  = P.T @  n_ax

        # re-write projected partials
        H_1_proj[k, 0:Npart] = np.hstack((h_p[0, :], r_p[0, :], n_p[0, 0]))
        H_2_proj[k, 0:Npart] = np.hstack((h_p[1, :], r_p[1, :], n_p[1, 1]))
        H_3_proj[k, 0:Npart] = np.hstack((h_p[2, :], r_p[2, :], n_p[2, 2]))

    print(" Finished!")

    # Write new regress with projected partials
    org_name = os.path.basename(reg_files[0])

    reg_out = outDir + '/' + org_name.replace('XX', 'NSM_1')
    print("Writing: " + os.path.basename(reg_out))
    reg_update_partials(reg_files[0],reg_out,np.arange(0,Npart - 2) + 1,H_1_proj)
    
    reg_out = outDir + '/' + org_name.replace('XX', 'NSM_2')
    print("Writing: " + os.path.basename(reg_out))
    reg_update_partials(reg_files[1],reg_out,np.arange(0,Npart - 2) + 1,H_2_proj)
    
    reg_out = outDir + '/' + org_name.replace('XX', 'NSM_3')
    print("Writing: " + os.path.basename(reg_out))
    reg_update_partials(reg_files[2],reg_out,np.arange(0,Npart - 2) + 1,H_3_proj)

    # end function
    return 0