from functions.read_reg_file_partials_batch import read_reg_batch
from functions.read_GG_obs import read_GG_obs
from functions.GG_rotation_partials import GG_rotation_partials
from functions.update_reg_partials_allRows_batch import reg_update_partials_batch

from tqdm import tqdm
import os
import shutil
import numpy as np

def modify_reg_NSM(GG_obs, reg_files, outDir, batch_size):
    """
    from the regress files for the XX, XY, XZ, YY, YZ, ZZ components and 
    the GG observations, construct the Null Space and modify the partials
    in the regress fields.

    Save the new regress files in the output folder location
    name: 

    Returns
    ----------
    N/A
    """

    # extract and load GG observations
    G          = read_GG_obs(GG_obs)

    Nt         = len(G)
    Nobs       = len(G.T)

    org_name     = os.path.basename(reg_files[0])
    reg_out_NSM1 = outDir + '/' + org_name.replace('XX', 'NSM1')
    reg_out_NSM2 = outDir + '/' + org_name.replace('XX', 'NSM2')
    reg_out_NSM3 = outDir + '/' + org_name.replace('XX', 'NSM3')

    # create a copy of the original regress files with the new name
    shutil.copy2(reg_files[0], reg_out_NSM1)
    shutil.copy2(reg_files[0], reg_out_NSM2)
    shutil.copy2(reg_files[0], reg_out_NSM3)

    # project partials
    print("Computing NSM projection to regress files")
    for k in tqdm(range(0, Nt, batch_size)):
        # Determine how many to read (handles the last smaller batch)
        current_batch = min(batch_size, Nt - k)

        # read original regress line by line
        H_XX = read_reg_batch(reg_files[0], start_obs=k, batch_size=batch_size)
        H_XY = read_reg_batch(reg_files[1], start_obs=k, batch_size=batch_size)
        H_XZ = read_reg_batch(reg_files[2], start_obs=k, batch_size=batch_size)
        H_YY = read_reg_batch(reg_files[3], start_obs=k, batch_size=batch_size)
        H_YZ = read_reg_batch(reg_files[4], start_obs=k, batch_size=batch_size)
        H_ZZ = read_reg_batch(reg_files[5], start_obs=k, batch_size=batch_size)

        Npart     = len(H_XX.T)

        # projected partials
        H_1_proj = np.zeros((batch_size, Npart))
        H_2_proj = np.zeros((batch_size, Npart))
        H_3_proj = np.zeros((batch_size, Npart))

        for i in range(current_batch):
            # partials
            h = np.array([H_XX[i,0:Npart-2],H_XY[i,0:Npart-2],H_XZ[i,0:Npart-2],\
                        H_YY[i,0:Npart-2],H_YZ[i,0:Npart-2],H_ZZ[i,0:Npart-2]])
            # residuals
            r = np.array([H_XX[i,Npart-2],H_XY[i,Npart-2],H_XZ[i,Npart-2],\
                        H_YY[i,Npart-2],H_YZ[i,Npart-2],H_ZZ[i,Npart-2] ]).reshape(6,1)
            # noise
            n = np.diag([H_XX[i,Npart-1], H_XY[i,Npart-1], H_XZ[i,Npart-1], \
                        H_YY[i,Npart-1], H_YZ[i,Npart-1], H_ZZ[i,Npart-1]])
            # observations
            g = G[k+i, 1:Nobs]

            # compute rotation partials
            Hrot = GG_rotation_partials(g)

            # SVD for the H_rot
            U, S, Vh = np.linalg.svd(Hrot.T)
            V = Vh.T # NOTE: Pyhton returns the Hermitian. We need to transpose

            # construct projector (Orthogonal version)
            P = V[:, 3:6]
            
            # project partials & residuals
            h_p  = P.T @  h
            r_p  = P.T @  r
            n_ax = n @ P
            n_p  = P.T @  n_ax

            # re-write projected partials 
            H_1_proj[i, 0:Npart] = np.hstack((h_p[0, :], r_p[0, :], n_p[0, 0]))
            H_2_proj[i, 0:Npart] = np.hstack((h_p[1, :], r_p[1, :], n_p[1, 1]))
            H_3_proj[i, 0:Npart] = np.hstack((h_p[2, :], r_p[2, :], n_p[2, 2]))

        # update regress files
        reg_update_partials_batch(
             reg_out_NSM1, 
             start_obs=k,  
             parCols=np.arange(0,Npart - 2)+1, 
             newVals_batch=H_1_proj,
             out_file=reg_out_NSM1
         )
        
        reg_update_partials_batch(
             reg_out_NSM2, 
             start_obs=k,  
             parCols=np.arange(0,Npart - 2)+1, 
             newVals_batch=H_2_proj,
             out_file=reg_out_NSM2
         )
        
        reg_update_partials_batch(
             reg_out_NSM3, 
             start_obs=k,  
             parCols=np.arange(0,Npart - 2)+1, 
             newVals_batch=H_3_proj,
             out_file=reg_out_NSM3
         )

    # end function
    return 0