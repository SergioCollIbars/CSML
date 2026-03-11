
from functions.read_reg_file_partials import read_reg
from functions.load_grav_field import load_matrix
from functions.read_apriori_reg import read_reg_apriori

import argparse
import numpy as np
import os

def main():
    # initialize parser
    parser = argparse.ArgumentParser()
    
    # 2. Add the arguments you want to receive
    # 'path' is the name of the variable, 'help' is what shows up if you type -h
    parser.add_argument("input_regress_folder", help="The path to the regress folder")
    args = parser.parse_args()

    file_folder = args.input_regress_folder
    INPUT_FOLDER = str(file_folder)
    all_files = os.listdir(INPUT_FOLDER)

    # script path
    script_dir         = os.path.dirname(os.path.abspath(__file__))

    # Point to the "data" folder
    data_folder        = os.path.join(script_dir, "data")
    GEO                = load_matrix(data_folder + "/GIF48.2007.txt")
    vals_GEO_zonal     = GEO[0:359, 0]
    val_GEO_sector     = GEO[359:len(GEO), :].reshape(-1, 1)

    vals_GEO = np.vstack((vals_GEO_zonal[0:4].reshape(4,1), val_GEO_sector))

    # read apriori values
    AP = read_reg_apriori(INPUT_FOLDER + "/" + all_files[0]) 

    # Information & Normal matrix
    Ax = 0; Nx = 0
    # read data
    for filename in all_files:
        # Check if it's a file and not a subfolder
        if os.path.isfile(os.path.join(INPUT_FOLDER, filename)):
            file        = INPUT_FOLDER + "/" + filename
            regres      = read_reg(file)

            Nparams     = len(regres.T)
            N           = Nparams -2 # Idex for max parameters in estimation. "Nparams-2" for all the SH partials

            residual    = regres[:, -2]
            sigma       = regres[:, -1]

            H           = regres[:, 0:N]
            Hb          = H[:, 0:6]                # bias partials
            Hpv         = H[:,6:12]                # pos & velocity partials
            Hcs         = H[:, 12:len(H)] * 1E6    # SH coefficient partials
            
            # update Normal equations
            W           = np.diag(1/(sigma**2))
            
            h           = np.hstack((Hb, Hcs))
            Ax          = Ax + h.T @ (W @ h)
            Nx          = Nx + h.T @ (W @ residual)

    # solve system of Linear Equations
    cov      = np.linalg.inv(Ax)
    cond_num = np.linalg.cond(cov)
    print(f"Condition number of the system: {cond_num}")

    # solve system
    x_hat    = cov @ Nx
    X_stm    = (AP[12:N] * 1E-6) + x_hat[6:len(x_hat)] 

    # read true coefficients
    SNR      = np.zeros((1, len(X_stm)))
    sig      = np.zeros((1, len(X_stm)))

    for i in range(0, len(X_stm)):
        e_c = vals_GEO[i] - X_stm[i]
        s_c = (cov[i][i])**(0.5)
        SNR[0,i]   = e_c / s_c
        sig[0,i]   = s_c

    print(f"Number of files processed: {len(all_files)}")
    
    # save computed data
    np.savetxt("data/SNR_grav_estim.txt", SNR, fmt='%.18e', delimiter='  ')
    np.savetxt("data/sigma_grav_estim.txt", sig, fmt='%.18e', delimiter='  ')

if __name__ == "__main__":
    main()