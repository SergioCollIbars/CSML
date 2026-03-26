import numpy as np

def GG_rotation_partials(g_obs):
    """
    Compute the rotation partials for a 3-2-1 rotation sequence
    (yaw - pitch -roll)

    Input
    -------
    g_obs: GG observations in the instrument frame. 
           Order: xx, xy, xz, yy, yz, zz

    Return
    -------
    H_rot: rotation partial matrix for yaw, pitch, roll errors
            Order: xx, xy, xz, yy, yz, zz
            Dimensions: 6 x 3
    """

    Yxx = g_obs[0]; Yxy = g_obs[1]; Yxz = g_obs[2]
    Yyy = g_obs[3]; Yyz = g_obs[4]; Yzz = g_obs[5]

    h = -np.array([[0, 2*Yxz, -2*Yxy], \
          [-Yxz, Yyz, Yxx - Yyy], \
          [Yxy, Yzz - Yxx, -Yyz], \
          [-2*Yyz, 0, 2*Yxy]    , \
          [Yyy - Yzz, -Yxy, Yxz], \
          [2*Yyz, -2*Yxz, 0]])
    
    # re-order to yaw, pitch, roll sequence
    # (h is in roll, pich, yaw sequence)
    Hrot = np.column_stack((h[:, 2], h[:, 1], h[:, 0]))

    return Hrot