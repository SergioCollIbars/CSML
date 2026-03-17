import numpy as np
from functions.read_att_files import read_att_to_DCM

def compute_att_error(file_AttObs_nom, file_AttObs_true, t0, tf):
    # read attitude files in Direction Cosine Matrices (DCM)
    DCM_nom  = read_att_to_DCM(str(file_AttObs_nom), t0, tf)
    DCM_true = read_att_to_DCM(str(file_AttObs_true), t0, tf)

    DCM_err = np.transpose(DCM_nom, (0, 2, 1)) @ DCM_true

    yaw   = np.arctan2(DCM_err[:, 0, 1], DCM_err[:, 0, 0])
    pitch = -np.arcsin(DCM_err[:, 2, 0])
    roll  = np.arctan2(DCM_err[:, 1, 2], DCM_err[:, 2, 2])

    # Stack into a (3, N) array to match your err_rpy shape
    err_rpy = np.array([roll, pitch, yaw])

    return err_rpy