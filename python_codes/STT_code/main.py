# This script is designed to propagate deviations to a spacecraft trajectory using higher order state transition tensors

# Arithmetic packages
import numpy as np

# Integration packages
from scipy.integrate import solve_ivp

# Plotting packages
import matplotlib.pyplot as plt

# Import other functions
from EOM import *
from STT_ODEs import *
from sim_param import *
from utility import *
from applySTTs import *

def main():
    # Compute the dx_k1 ... dx_dp terms in the solution flow expansion
    dx_list = compute_dx(dX0, order)

    # Integrate spacecraft trajectory and state transition tensors evaluated along the trajectory
    solrk = solve_ivp(STT_ODEs, [t[0],t[-1]], IC, 'DOP853', t, args = (param,dynamics), rtol=1e-13, atol=1e-13)

    # Extract state and STT information to propagate trajectory
    dX_hist = np.zeros((len(t),stateSize))
    X_hist = np.zeros((len(t),stateSize))
    for ii in range(len(t)):
        traji =  solrk.y[:,ii]
        Xi, phi_i_list = reshape_STT(traji,order,stateSize) # phi_i_list contains the STM, 2nd order STT, the 3rd order STT, ...
        dXi = STT_update(phi_i_list, dx_list, order, stateSize)
        dX_hist[ii,:] = dXi
        X_hist[ii,:] = Xi

    # Reconstruct the perturbed trajectory using the mapped deviations
    Xpert_hist = X_hist + dX_hist

    # Plot perturbed trajectory approximated with STTs, perturbed trajectory with the full nonlinear ODEs, and their error
    solrk_pert = solve_ivp(dynamics, [t[0],t[-1]], X0 + dX0, 'DOP853', t, args = (param,), rtol=1e-13, atol=1e-13)

    fig, (ax1, ax2, ax3) = plt.subplots(1, 3)
    plot1 = ax1.scatter(Xpert_hist[:,0], Xpert_hist[:,1])
    plot2 = ax2.scatter(solrk_pert.y[0], solrk_pert.y[1])
    plot3 = ax3.scatter(t, np.linalg.norm(Xpert_hist - solrk_pert.y.T, axis = 1))
    ax1.set_xlabel('X [nd]')
    ax1.set_ylabel('Y [nd]')
    ax1.set_title(sim_type+': Order '+str(order)+' Propagation')
    ax1.grid()
    ax2.set_xlabel('X [nd]')
    ax2.set_ylabel('Y [nd]')
    ax2.set_title(sim_type+': Full Nonlinear Propagation')
    ax2.grid()
    ax3.set_xlabel('Time [nd]')
    ax3.set_ylabel(r'$|X_{STT} - X_{NL}|$ [nd]')
    ax3.set_title(sim_type+': Error')
    ax3.set_yscale('log')
    ax3.grid()    
    plt.show()
    plt.close()

main()