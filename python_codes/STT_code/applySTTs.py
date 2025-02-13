# This scripts contains functions designed to compute the terms in the higher order solution flow expansion

# Arithmetic packages
import numpy as np

# ----------------------------------------------------------------------------
# The following function computes the state deviations up to 4th order; 
# That is, this function computes the dx_k1 ... dx_kp terms
# ----------------------------------------------------------------------------
def compute_dx(dx, order):
    # dx -> spacecraft state perturbation
    dx_list = []

    if order >= 1:
        # FIRST ORDER
        # --------------------------------------------------------
        dx_first = dx
        dx_list.append(dx_first)

    if order >= 2:      
        # SECOND ORDER
        # --------------------------------------------------------
        dx_second = np.einsum('i,j->ij',dx,dx)
        dx_list.append(dx_second)
            
    if order >= 3:
        # THIRD ORDER
        # --------------------------------------------------------
        dx_third = np.einsum('i,j,k->ijk',dx,dx,dx,optimize=True)
        dx_list.append(dx_third)
    
    if order >= 4:    
        # FOURTH ORDER
        # --------------------------------------------------------
        dx_fourth = np.einsum('i,j,k,l->ijkl',dx,dx,dx,dx,optimize=True)

        dx_list.append(dx_fourth)
    return dx_list

# ----------------------------------------------------------------------------
# The following function computes the update to the state deviation using
# higher order state transition tensors
# ----------------------------------------------------------------------------
def STT_update(phi_list, dx_list, order, stateSize):
    dx = np.zeros(stateSize)

    if order >= 1:
        dx += np.einsum('ia,a->i',phi_list[0],dx_list[0])
        
        if order >= 2:
            dx += 1/2*np.einsum('iab,ab->i',phi_list[1],dx_list[1],optimize=True)
            
            if order >= 3:
                dx += 1/6*np.einsum('iabc,abc->i',phi_list[2],dx_list[2],optimize=True)
                
                if order >= 4:
                    dx += 1/24*np.einsum('iabcd,abcd->i',phi_list[3],dx_list[3],optimize=True)

    return dx