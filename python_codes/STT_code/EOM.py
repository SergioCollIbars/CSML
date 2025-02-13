# This script will contain a collection of all the dynamical systems we wish to consider for our spacecraft

# Arithmetic packages
import numpy as np

# Dynamical Models
def CR3BP(t,X,param):
    # System Mass Parameters
    mu = param[0]

    # Unpack the spacecraft state
    x = X[0]
    y = X[1]
    z = X[2]
    xdot = X[3]
    ydot = X[4]
    zdot = X[5]
    
    rho1 = (y**2 + z**2 + (-1 + x + mu)**2)
    rho2 = (y**2 + z**2 + (x + mu)**2)
    
    # Compute accelerations
    xdotdot = 2*ydot + x - (mu*(-1 + x + mu))/rho1**(3/2) - (1-mu)*(x+mu)/rho2**(3/2)
    ydotdot = -2*xdot + y - y*mu/rho1**(3/2) - y*(1-mu)/rho2**(3/2)
    zdotdot = -z*mu/rho1**(3/2) - z*(1-mu)/rho2**(3/2)
    
    
    dX = np.array([xdot, ydot, zdot, xdotdot, ydotdot, zdotdot])
    return dX