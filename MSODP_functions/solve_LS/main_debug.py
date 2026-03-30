
from functions.read_reg_file_partials import read_reg
from functions.load_grav_field import load_matrix
from functions.read_apriori_reg import read_reg_apriori

import numpy as np
import os

# Folder with regress files
INPUT_FOLDER = "/Users/sergiocollibars/Documents/regres_2008_5by5"

def main():
    all_files = [f for f in os.listdir(INPUT_FOLDER) if os.path.isfile(os.path.join(INPUT_FOLDER, f))]

    # 2. Load Reference/GEO data
    # (Assuming these helper functions are defined elsewhere in your script)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    data_folder = os.path.join(script_dir, "data")
    
    # Load GEO matrix and slice values
    GEO = load_matrix(os.path.join(data_folder, "GIF48.2007.txt"))
    vals_GEO_zonal = GEO[0:359, 0]
    val_GEO_sector = GEO[359:len(GEO), :].reshape(-1, 1)
    vals_GEO = np.vstack((vals_GEO_zonal[0:4].reshape(4,1), val_GEO_sector))

    # Read apriori values from the first file
    AP = read_reg_apriori(os.path.join(INPUT_FOLDER, all_files[0])) 

    # 3. First Pass: Accumulate Normal Equations & Prefit
    prefit = [] 
    Ax = None
    Nx = None
    valid_files = [] # Store files that were successfully read

    # sigma scale value
    # sigma_scale = 3E-4
    sigma_scale = 1

    print("Processing files for Normal Equations...")
    for filename in all_files:
        file_path = os.path.join(INPUT_FOLDER, filename)
        regres = read_reg(file_path)
        
        # Slicing logic
        Nparams = regres.shape[1]
        N = Nparams - 2 
        
        residual = regres[:, -2]
        sigma = regres[:, -1] * sigma_scale # Weighting factor

        print(f'Mean sigma {np.max(sigma)}')
        
        # Design Matrix H
        H   = regres[:, 0:N]
        Hb  = H[:, 0:6]                # bias partials
        Hcs = H[:, 12:] * 1E6          # SH coefficient partials
        h   = np.hstack((Hb, Hcs))
        
        # Memory-efficient weighting
        W_vec = 1 / (sigma**2)
        
        # Accumulate Ax and Nx
        # (h.T * W_vec) @ h is equivalent to h.T @ diag(W) @ h but faster
        current_Ax = (h.T * W_vec) @ h
        current_Nx = (h.T * W_vec) @ residual

        if Ax is None:
            Ax = np.zeros_like(current_Ax)
            Nx = np.zeros_like(current_Nx)

        Ax += current_Ax
        Nx += current_Nx
        
        # Save prefit residuals
        prefit.extend(residual.tolist())
        valid_files.append(file_path)

    # 4. Solve the System
    # Ax * x_hat = Nx  => x_hat = inv(Ax) * Nx
    cov = np.linalg.inv(Ax)
    x_hat = cov @ Nx
    
    # Condition number check
    cond_num = np.linalg.cond(Ax) 
    print(f"Condition number of the system: {cond_num:.2e}")

    # 5. Second Pass: Compute Postfit Residuals
    postfit = []
    print("Computing postfit residuals...")
    for file_path in valid_files:
        regres = read_reg(file_path)
        residual = regres[:, -2]
        
        H = regres[:, 0:(regres.shape[1]-2)]
        h = np.hstack((H[:, 0:6], H[:, 12:] * 1E6))
        
        # Postfit = Observation - Prediction
        # current_postfit = residual - (H * delta_x)
        current_postfit = residual - (h @ x_hat)
        postfit.extend(current_postfit.tolist())

    # 6. Calculate Statistics
    prefit_array = np.array(prefit)
    postfit_array = np.array(postfit)

    rms_pre = np.sqrt(np.mean(prefit_array**2))
    rms_post = np.sqrt(np.mean(postfit_array**2))

    print("-" * 30)
    print(f"RMS prefit value:  {rms_pre:.6e} 1/s^2")
    print(f"RMS postfit value: {rms_post:.6e} 1/s^2")
    print(f"Improvement:       {((rms_pre - rms_post)/rms_pre)*100:.2f}%")
    print("-" * 30)

    # 7. SNR and Coefficient Estimation
    X_stm = (AP[12:12+len(x_hat)-6] * 1E-6) + x_hat[6:] 

    SNR = np.zeros((1, len(X_stm)))
    sig = np.zeros((1, len(X_stm)))
    err = np.zeros((1, len(X_stm)))

    for i in range(len(X_stm)):
        e_c = vals_GEO[i] - X_stm[i]
        s_c = np.sqrt(cov[i+6, i+6]) # Variance from covariance matrix diagonal
        SNR[0, i] = e_c[0] / s_c
        sig[0, i] = s_c
        err[0, i] = e_c[0]

    # 8. Save output
    if not os.path.exists("data"):
        os.makedirs("data")
        
    np.savetxt("data/SNR_grav_estim.txt", SNR, fmt='%.18e', delimiter='  ')
    np.savetxt("data/sigma_grav_estim.txt", sig, fmt='%.18e', delimiter='  ')
    np.savetxt("data/err_grav_estim.txt", err, fmt='%.18e', delimiter='  ')
    
    print(f"Processing complete. Files processed: {len(valid_files)}")

if __name__ == "__main__":
    main()