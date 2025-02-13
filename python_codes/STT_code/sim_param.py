# This script will contain the simulation parameters to model the spacecraft evolution

# Arithmetic packages
import numpy as np

# Import the equations of motion
from EOM import *

sim_type = "CR3BP"
if sim_type == "CR3BP":
    # Earth-Moon CR3BP
    dynamics = CR3BP
    GM_E = 398600.435436 #km^3/s^2
    GM_M = 4902.800066 #km^3/s^2
    mu = GM_M/(GM_M + GM_E)

    # Dimensional scaling units
    l_EM = 384400 #in km; Average Earth-Moon separation i.e semi-major axis
    t_EM = np.sqrt(l_EM**3/(GM_E+GM_M)) #non-dimensional time unit in seconds

    # L2 Lyapunov Orbit from JPL Three-Body Periodic Orbits Library
    # x0 = 1.0953533743235189E+0
    # y0 = -1.0879975950267760E-28
    # z0 = 0
    # x_dot0 = 1.3016066486537214E-15
    # y_dot0 = 2.9531900698678965E-1
    # z_dot0 = 0
    # t0 = 0
    # tf = 3.4981205539386764E+0
    # N = 1000
    # t = np.linspace(t0,tf,N)
    # X0 = np.array([x0, y0, z0, x_dot0, y_dot0, z_dot0])

    # L1 NHRO orbit
    x0 = 1.021968177072928
    y0 = 0
    z0 = -0.18206
    x_dot0 = 0
    y_dot0 = -0.1031401430288178
    z_dot0 = 0
    t0 = 0
    tf = 1.3819
    N = 1000
    t = np.linspace(t0,tf,N)
    X0 = np.array([x0, y0, z0, x_dot0, y_dot0, z_dot0])

    # Initial state perturbation
    dX0 = np.array([1E-6,0,0,0,0,0])
    
    # STT parameters
    order = 2 # Truncate after 2nd order expansion
    stateSize = 6

    # Simulation parameters
    param = np.array([mu, order, stateSize])

    # Initialize state and higher order expansions
    phi_0_first = np.eye(stateSize)
    phi_0_second = np.zeros((stateSize,stateSize,stateSize))
    phi_0_third = np.zeros((stateSize,stateSize,stateSize,stateSize))
    phi_0_fourth = np.zeros((stateSize,stateSize,stateSize,stateSize,stateSize))

    # Formulate full initial state vector
    IC = np.hstack([X0, phi_0_first.reshape(-1)])

    if order >= 2:
        IC = np.hstack([IC, phi_0_second.reshape(-1)])
        
        if order >= 3:
            IC = np.hstack([IC, phi_0_third.reshape(-1)])
            
            if order >= 4:
                IC = np.hstack([IC, phi_0_fourth.reshape(-1)])