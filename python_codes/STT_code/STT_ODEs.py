# This script will compute the ordinary differential equations for state transition tensors for a given dynamical system

# Arimethic packages
import numpy as np

# Automatic differentiation packages
from pyaudi import gdual_double as gdual
from pyaudi import sqrt

def STT_ODEs(t, X, param, dynamics):
    # Automatic differentiation parameters
    order = int(param[1]) # MUST be an integer
    stateSize = int(param[2]) # MUST be an integer

    # Reshape STM from 36 x 1 vector to 6 x 6 matrix, 2nd order STT from 216 x 1 vector 6 x 6 x 6 rank 3 tensor, ...
    if order >= 1:
        second_start = int(stateSize + stateSize**2)
        phi_first = X[stateSize:second_start].reshape((stateSize,stateSize))
        
        third_start = second_start + stateSize**3
        if order >= 2:
            phi_second = X[second_start:third_start].reshape((stateSize,stateSize,stateSize))
            
            if order >= 3:
                fourth_start = third_start + stateSize**4
                phi_third = X[third_start:fourth_start].reshape((stateSize,stateSize,stateSize,stateSize))
                
                if order >= 4:
            
                    phi_fourth = X[fourth_start:].reshape((stateSize,stateSize,stateSize,stateSize,stateSize))
    
    
    # Compute A matrices/tensors (partial derivatives of our dynamics)
    if order >= 1:
        A_first = np.zeros((stateSize,stateSize))
        A_second = np.zeros((stateSize,stateSize,stateSize))
        A_third = np.zeros((stateSize,stateSize,stateSize,stateSize))
        A_fourth = np.zeros((stateSize,stateSize,stateSize,stateSize,stateSize))
        
        # Initialize state as gdual variables
        g_order = order
        X_i_g = [gdual(X[0],"x", g_order), gdual(X[1], "y", g_order), gdual(X[2], "z", g_order), 
                 gdual(X[3],"xdot", g_order), gdual(X[4], "ydot", g_order), gdual(X[5], "zdot", g_order)]
        
        # Compute dynamics with gdual variables
        F = dynamics(t,X_i_g, param)
        dual_str = ["x","y","z","xdot","ydot","zdot"]
        
        
        # Extract partial derivatives using automatic differentiation package
        for i in range(stateSize):
            for a in range(stateSize):
               A_first_gdual = F[i].partial(dual_str[a])
               A_first[i,a] = A_first_gdual.constant_cf
               
               
               if order >= 2:
                   for b in range(stateSize):
                        A_second_gdual = A_first_gdual.partial(dual_str[b])
                        A_second[i,a,b] = A_second_gdual.constant_cf
                        
                        
                        if order >= 3:
                            for c in range(stateSize):
                                
                                A_third_gdual = A_second_gdual.partial(dual_str[c])
                                A_third[i,a,b,c] = A_third_gdual.constant_cf
                               
                                if order >= 4:
                                    for d in range(stateSize):
                                        A_fourth_gdual = A_third_gdual.partial(dual_str[d])
                                        A_fourth[i,a,b,c,d] = A_fourth_gdual.constant_cf
                                
    
    # -------------------------------------------------------------------------                
    
    # Compute time derivatives of the STM, 2nd order STT, ...
    # First order
    if order >= 1:
        phidot_first = np.einsum('ij,ja->ia',A_first,phi_first)
        phidot = phidot_first.reshape(-1)
    else:
        phidot = []

    # Second order
    if order >= 2:
        phidot_second = np.zeros((stateSize,stateSize,stateSize))
                          
        phidot_second += np.einsum('ij,jab->iab',A_first,phi_second)
        phidot_second += np.einsum('ijk,ja,kb->iab',A_second,phi_first,phi_first,optimize=True)
        
        phidot = np.hstack([phidot, phidot_second.reshape(-1)])
     
                                       
    # Third order
    if order >= 3:
        phidot_third = np.zeros((stateSize,stateSize,stateSize,stateSize)) 
                      
        phidot_third += np.einsum('ij,jabc->iabc',A_first,phi_third)
        phidot_third += np.einsum('ijk,ja,kbc->iabc',A_second,phi_first,phi_second,optimize=True) 
        phidot_third += np.einsum('ijk,jab,kc->iabc',A_second,phi_second,phi_first,optimize=True)
        phidot_third += np.einsum('ijk,jac,kb->iabc',A_second,phi_second,phi_first,optimize=True)
        phidot_third += np.einsum('ijkl,ja,kb,lc->iabc',A_third,phi_first,phi_first,phi_first,optimize=True)
        
        phidot = np.hstack([phidot, phidot_third.reshape(-1)])
    

    # Fourth order   
    if order >= 4:
        phidot_fourth = np.zeros((stateSize,stateSize,stateSize,stateSize,stateSize))
        
        phidot_fourth += np.einsum('ij,jabcd->iabcd',A_first,phi_fourth)
        phidot_fourth += np.einsum('ijk,jabc,kd->iabcd',A_second,phi_third,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijk,jabd,kc->iabcd',A_second,phi_third,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijk,jacd,kb->iabcd',A_second,phi_third,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijk,jab,kcd->iabcd',A_second,phi_second,phi_second,optimize=True)
        phidot_fourth += np.einsum('ijk,jac,kbd->iabcd',A_second,phi_second,phi_second,optimize=True)
        phidot_fourth += np.einsum('ijk,jad,kbc->iabcd',A_second,phi_second,phi_second,optimize=True)
        phidot_fourth += np.einsum('ijk,ja,kbcd->iabcd',A_second,phi_first,phi_third,optimize=True)
        phidot_fourth += np.einsum('ijkl,jab,kc,ld->iabcd',A_third,phi_second,phi_first,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijkl,jac,kb,ld->iabcd',A_third,phi_second,phi_first,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijkl,jad,kb,lc->iabcd',A_third,phi_second,phi_first,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijkl,ja,kbc,ld->iabcd',A_third,phi_first,phi_second,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijkl,ja,kbd,lc->iabcd',A_third,phi_first,phi_second,phi_first,optimize=True)
        phidot_fourth += np.einsum('ijkl,ja,kb,lcd->iabcd',A_third,phi_first,phi_first,phi_second,optimize=True)
        phidot_fourth += np.einsum('ijklm,ja,kb,lc,md->iabcd',A_fourth,phi_first,phi_first,phi_first,phi_first,optimize=True)
        
        phidot = np.hstack([phidot, phidot_fourth.reshape(-1)])
        
    
    # Reshape velocities, accelerations, and STT time derivatives into a column vector
    dX = np.hstack([dynamics(t,X[0:stateSize], param),phidot])
    return dX