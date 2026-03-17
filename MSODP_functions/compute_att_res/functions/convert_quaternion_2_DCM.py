import numpy as np

def quat2dcm(q):
    """
    Translates a single [q0, q1, q2, q3] quaternion to a 3x3 DCM.
    """
    q = np.asarray(q).flatten()
    if len(q) != 4:
        raise ValueError("Input quaternion must have 4 elements.")
    
    # Normalize
    q = q / np.linalg.norm(q)
    q0, q1, q2, q3 = q

    # Construct the matrix
    C = np.array([
        [1 - 2*(q2**2 + q3**2), 2*(q1*q2 - q0*q3),     2*(q1*q3 + q0*q2)],
        [2*(q1*q2 + q0*q3),     1 - 2*(q1**2 + q3**2), 2*(q2*q3 - q0*q1)],
        [2*(q1*q3 - q0*q2),     2*(q2*q3 + q0*q1),     1 - 2*(q1**2 + q2**2)]
    ])
    return C