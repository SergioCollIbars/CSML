from functions.modify_reg_NSM import modify_reg_NSM
from functions.read_reg_file_partials import read_reg

# Input / output files
inReg  = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/regress_reg/goce_XX_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg"
outReg = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/NSM_regress_reg"

inObs  = "/Users/sergiocollibars/Desktop/CSML/MSODP_functions/GG_apriori/goce_eggreg_2008-08-01_RL5061_RL05surpv6.980.001.ggr"

reg_files = [inReg, inReg, inReg, inReg, inReg, inReg]

# run NSM transformation in the regress files
Nbatch_size = 3
modify_reg_NSM(inObs, reg_files, outReg, Nbatch_size)

# read one of the output generated files
H_NSM1 = read_reg(outReg + "/goce_NSM1_eggreg_2008-08-01_RL5061_RL05_v6.980_10uE.001.reg")

# display first N rows and first M columsn
Nrows = 10; Ncols = 6
print(H_NSM1[0:Nrows,0:Ncols])