from functions.read_reg_file_partials import read_reg
from functions.update_reg_partials_allRows import reg_update_partials

import numpy as np
import subprocess

"""
TEST READ AND WRITE REGRESS ROUTINES

Author: Sergio Coll-Ibars
Date:   05/03/2026

Description: Test the routines used in v0 code version to
    write and write the regress files.
    Use the 'dmpreg.c' code developed at CSR to convert 
    regres .reg into ASCII format. 
    Compare outputs and check the differences 
"""

INPUT_REG  = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/regres_files/" +\
             "goce_XX_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg"
OUTPUT_REG = "./data/goce_NSM_1_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg"

OUTPUT_ASCII_ORG  = "./data/original.out"
OUTPUT_ASCII_TEST = "./data/test.out"

# read the partials values in the regress files
H  = read_reg(INPUT_REG)

parCols = np.arange(0, (len(H.T) - 2)) + 1 
Nobs    = len(H)

Hc = np.zeros((Nobs ,len(H.T)))

# check the Identity matrix
Nbatch  = 6
I  = np.eye(Nbatch)

for k in range(0, Nobs, Nbatch):
    batch_end = min(k + Nbatch, Nobs)
    val       = H[k:batch_end, :]
    val_new   = I @ val
    Hc[k:batch_end, :] = val_new 

# make sure the XX, XY, XZ .. are well asigned to the variables
# check the ordering of matrices (Row-Column indexing)

# # Modify only some values
# Hc[2, 0] = 0

# write the partials values back
reg_update_partials(INPUT_REG, OUTPUT_REG, parCols, Hc)

# use dmpreg to generate ASCII files
code   = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/CSR_lonestar6_codes/dmpreg_code/dmpreg"
command = f"{code} {INPUT_REG} 3 -1 > {OUTPUT_ASCII_ORG}"
subprocess.run(command, shell=True, executable="/bin/zsh", check=True)

command = f"{code} {OUTPUT_REG} 3 -1 > {OUTPUT_ASCII_TEST}"
subprocess.run(command, shell=True, executable="/bin/zsh", check=True)

# check differnces between both and print
char_limit = 100
command = f"diff {OUTPUT_ASCII_ORG} {OUTPUT_ASCII_TEST} | cut -c 1-{char_limit}"
subprocess.run(command, shell=True)