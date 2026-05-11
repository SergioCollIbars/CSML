from functions.read_reg_file_partials import read_reg
from functions.update_reg_partials_allRows import reg_update_partials
import numpy as np
from pathlib import Path
import os

"""
MODIFY REGRES FILES MEAS. WEIGTH

Description: modify the measurement weigth in the regres files to all 
            the files in the specified folder
"""

def main():
    folder_reg   = Path("/Users/sergiocollibars/Documents/regres_files").resolve()
    meas_weigth  = 1E-6
    folder = Path("output")
    folder.mkdir(parents=True, exist_ok=True)

    # Iterate through OBS files
    for reg_file in sorted(folder_reg.glob("*.reg")):  # Adjust extension if needed

        H = read_reg(reg_file)
        Npart     = len(H.T)

        # get measurement weigth (normalized)
        weigth = H[:, -1]
        weigth_norm = weigth / weigth

        # new measurement weigth
        weigth_new  = weigth_norm * meas_weigth

        # re-write measurement weigth
        H[:, -1] = weigth_new
        org_name = os.path.basename(reg_file)
        reg_file_out = "output/" + org_name
        reg_update_partials(reg_file, reg_file_out ,np.arange(0,Npart - 2) + 1, H)
        
if __name__ == "__main__":
    main()