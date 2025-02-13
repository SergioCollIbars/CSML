#---------------------------------------------------------------------#
#                       INITIAL APPROACH TO QGG                       # 
#                                                                     #
# Author: Sergio Coll                                                 #
# Date: 08/25/22                                                      #
# Description: Solve the problem of a circular orbit around Earth     #
#               Get the measure obtained by the QGG in this orbit     #
#---------------------------------------------------------------------#


# --------------------------------------------------------------------#
#                               IMPORTS                               #
# --------------------------------------------------------------------#
from cmath import pi
import constants as c
import functions as f
import numpy as np
import sys
import math 
from tqdm import tqdm

# --------------------------------------------------------------------#
#                            CONFIGURATION                            #
# --------------------------------------------------------------------#

mainConfig = int(sys.argv[1])               # Main configuration [1 / 2]

gravityPotential = 1                        # Number of terms used to express potential gravity.
latitude = 90                               # Orbit latitude [deg]
longitude = 90                              # Orbit longitude [deg]
r = 1000e3                                  # Orbit height over planet surface [m]
A1 = np.array([[0.5, 0 , 0]])               # QGG distance separation x, y, z axis along COM sensor 1 [m]
A2 = np.array([[-0.5, 0 , 0]])              # QGG distance separation x, y, z axis along COM sensor 2 [m]
gamma0 = 0                                  # Initial swipe angle [rad]
gamma1 = pi / 2                             # Final swipe angle [rad] 

points = 100                                # Number of point in plots

# --------------------------------------------------------------------#
#                               MAIN1()                               #
# --------------------------------------------------------------------#

def main1():

    # Decription: Use matrix equation formulation to asses QGG mesuraments at different
    #   positions

    # Obtain gravity values:
    r_o = np.array([[c.R_E + r], [0],[0]])

    # U = f.getPotential(gravityPotential, c.G, c.M_E, r_o, longitude, latitude)
    [Omega_ox, Omega_oy, Omega_oz] = f.getAngularVel(c.G, c.M_E, f.norm(r_o))
    omega = np.array([[0], [Omega_ox], [Omega_oy], [Omega_oz]])

    M  = np.array([[0], [0], [0], [0]])                                                 # Input torque to the satellite
    omegaDot = np.array(np.zeros((4, 1)))
    omegaDot= f.getAngularAcc(omega, omegaDot, M)
    OmegaDot_ox = omegaDot[1][0]
    OmegaDot_oy = omegaDot[2][0]
    OmegaDot_oz = omegaDot[3][0]

    # Construct angular acceleration matrixs
    Omega2 = np.array([[-pow(Omega_oz,2) -pow(Omega_oy,2), Omega_ox * Omega_oy, Omega_ox * Omega_oz],
              [Omega_ox * Omega_oy, -pow(Omega_oz,2) -pow(Omega_ox,2), Omega_oy * Omega_oz],
              [Omega_ox * Omega_oz, Omega_oy * Omega_oz, -pow(Omega_ox,2) -pow(Omega_oy,2)]])
    
    OmegaDot = np.array([[0, -OmegaDot_oz, OmegaDot_oy],
                [OmegaDot_oz, 0, -OmegaDot_ox],
                [-OmegaDot_oy, OmegaDot_ox, 0]])
    
    
    # Compute QGG acceleration (difference between two accelerometers)
    angle = np.linspace(gamma0 ,gamma1, num=points)
    ad = np.zeros([3, len(angle)])
    cont = 0
    for k in angle:
        Delta_r1 = f.norm(A1.T) * np.array([[np.cos(k)],                        # Sensor 1 radius, Axis: {xs, ys, zs}
                                [np.sin(k)],
                                [0]])
        Delta_r2 = -f.norm(A2.T) * np.array([[np.cos(k)],                       # Sensor 2 radius, Axis: {xs, ys, zs}
                                [np.sin(k)],
                                [0]])
            

        U = np.array([[2 * c.G * c.M_E / pow(f.norm(r_o), 3), 0, 0],
         [0, -c.G * c.M_E / pow(f.norm(r_o), 3), 0],
         [0, 0, -c.G * c.M_E / pow(f.norm(r_o), 3)]])

        # Compute acceleration
        a1 = np.dot(-(U - Omega2 - OmegaDot), Delta_r1)
        a2 = np.dot(-(U - Omega2 - OmegaDot), Delta_r2) 

        ad[0][cont] = a1[0] - a2[0]
        ad[1][cont] = a1[1] - a2[1]
        ad[2][cont] = a1[2] - a2[2]

        cont = cont + 1

    # Plot results
    f.plotAcc(angle, ad, '/Users/sergiocollibars/Desktop/CSML/figures/main1_fig')
    
    
    return 0


# --------------------------------------------------------------------#
#                               MAIN2()                               #
# --------------------------------------------------------------------#

def main2():

    # Decription: Use vector equation formulation to asses QGG mesuraments at different
    #   positions
    
    # Compute radius 
    r_o = np.array([[c.R_E + r],                                                 # Orbital radius, Axis: {xi, yi, zi}
                    [0],
                    [0]])
    Delta_r1 = f.norm(A1.T) * np.array([[np.cos(gamma0)],                        # Sensor 1 radius, Axis: {xs, ys, zs}
                                [np.sin(gamma0)],
                                [0]])
    Delta_r2 = -f.norm(A2.T) * np.array([[np.cos(gamma0)],                       # Sensor 2 radius, Axis: {xs, ys, zs}
                                [np.sin(gamma0)],
                                [0]])
    # Compute orbital values
    [Omega_ox, Omega_oy, Omega_oz] = f.getAngularVel(c.G, c.M_E, f.norm(r_o))
    Omega = np.array([[Omega_ox],                                               # Angular velocity, Axis: {xs, ys, zs}
                      [Omega_oy],
                      [Omega_oz]])
    
    mu = c.G * c.M_E                                                            # Mu parameter

    # Compute accelerometer meassurements: a1, a2
    # Guess acceleration at differnent sensor positions

    angle = np.linspace(gamma0 ,gamma1, num=points)
    ad = np.zeros([3, len(angle)])
    cont = 0
    for i in angle:
        Delta_r1 = f.norm(A1.T) * np.array([[np.cos(i)],                       # Sensor 1 radius, Axis: {xs, ys, zs}
                                [np.sin(i)],
                                [0]])
        Delta_r2 = -f.norm(A2.T) * np.array([[np.cos(i)],                      # Sensor 2 radius, Axis: {xs, ys, zs}
                                [np.sin(i)],
                                [0]])

        a1 = np.cross(Omega.T, np.cross(Omega.T, r_o.T)).T + np.cross(Omega.T, np.cross(Omega.T, Delta_r1.T)).T + mu / pow(f.norm(r_o + Delta_r1), 3) * (r_o + Delta_r1)
        a2 = np.cross(Omega.T, np.cross(Omega.T, r_o.T)).T + np.cross(Omega.T, np.cross(Omega.T, Delta_r2.T)).T + mu / pow(f.norm(r_o + Delta_r2), 3) * (r_o + Delta_r2)
        
       #  print(f'iteration {cont} \n a2:{a2} \n a1: {a1}')

        ad[0][cont] = a1[0] - a2[0]
        ad[1][cont] = a1[1] - a2[1]
        ad[2][cont] = a1[2] - a2[2]
        cont = cont + 1
    
    # Plot results
    print(f'Theorical max acc val: {-6 * mu / pow(f.norm(r_o),3) * f.norm(A1.T)}')
    f.plotAcc(angle, ad, '/Users/sergiocollibars/Desktop/CSML/figures/main2_fig')
    

    return 0


# --------------------------------------------------------------------#
#                               MAIN3()                               #
# --------------------------------------------------------------------#

# WARNING: Needs to be re-coded and adapted to the new structure!
        
def main3():

    # Description: Using the Matlab SS code (altitude determination code) evaluate how the QGG
    #   mesuraments are afected. 

    # Obtain gravity values:
    r_o = np.array([[c.R_E + r], [0],[0]])

    # Read data from txt files
    omega = f.readFile('/Users/sergiocollibars/Desktop/CSML/codes/matlab_datafiles/angular_vel_results.txt')
    omegaDot = f.readFile('/Users/sergiocollibars/Desktop/CSML/codes/matlab_datafiles/angular_acc_results.txt')
    moments = f.readFile('/Users/sergiocollibars/Desktop/CSML/codes/matlab_datafiles/moments_results.txt')
    time = moments[:][0]

    def computeAcc(Omega_ox, Omega_oy, Omega_oz, OmegaDot_ox, OmegaDot_oy, OmegaDot_oz):

        # Compute gravity gradient tensor
        U = np.array([[2 * c.G * c.M_E / pow(f.norm(r_o), 3), 0, 0],
            [0, -c.G * c.M_E / pow(f.norm(r_o), 3), 0],
            [0, 0, -c.G * c.M_E / pow(f.norm(r_o), 3)]])

        # Construct angular acceleration matrixs
        Omega2 = np.array([[-pow(Omega_oz,2) -pow(Omega_oy,2), Omega_ox * Omega_oy, Omega_ox * Omega_oz],
                [Omega_ox * Omega_oy, -pow(Omega_oz,2) -pow(Omega_ox,2), Omega_oy * Omega_oz],
                [Omega_ox * Omega_oz, Omega_oy * Omega_oz, -pow(Omega_ox,2) -pow(Omega_oy,2)]])
    
        OmegaDot = np.array([[0, -OmegaDot_oz, OmegaDot_oy],
                    [OmegaDot_oz, 0, -OmegaDot_ox],
                    [-OmegaDot_oy, OmegaDot_ox, 0]])

        # Compute distance
        Delta_r1 = f.norm(A1.T) * np.array([[np.cos(gamma0)],                        # Sensor 1 radius, Axis: {xs, ys, zs}
                                [np.sin(gamma0)],
                                [0]])
        Delta_r2 = -f.norm(A2.T) * np.array([[np.cos(gamma0)],                       # Sensor 2 radius, Axis: {xs, ys, zs}
                                [np.sin(gamma0)],
                                [0]])

        # Compute acceleration
        a1 = np.dot(-(U - Omega2 - OmegaDot), Delta_r1)
        a2 = np.dot(-(U - Omega2 - OmegaDot), Delta_r2) 
        return (a1 - a2) # x, y, z

    ad = np.array(np.zeros((4, 3,len(omega))))
    for t in range(len(omega)):

        p0 =  math.sqrt(c.G * c.M_E / pow(r, 3))                                        # Orbital angular velocity   [rad / s]

        # compute angular spacecraft displacement at current time
        gamma = omega[t][3] * time[t]
        beta =  omega[t][2] * time[t]
        theta = omega[t][1] * time[t]

        # Compute rotational matrix {x', y', z'} -> {xs, ys, zs}
        R_z = np.array([[np.cos(gamma), -np.sin(gamma), 0],                        
                        [np.sin(gamma), np.cos(gamma), 0],
                        [0, 0, 1]])
        R_y = np.array([[np.cos(beta), 0, np.sin(beta)],
                        [0, 1, 0],
                        [-np.sin(beta), 0, np.cos(beta)]])
        R_x = np.array([[1, 0, 0],
                        [0, np.cos(theta), -np.sin(theta)],
                        [0, np.sin(theta), np.cos(theta)]])

        R = np.dot(R_x, np.dot(R_y, R_z))

        r_i = f.norm(r_o) * np.array([[np.cos(p0 * time[t])], [-np.sin(p0 * time[t])], [0]])        # inertial frame radious
        r_b = R.T * r_i                                                                             # body frame radious

        [OmegaDot_ox, OmegaDot_oy, OmegaDot_oz] = f.getAngularAcc(0, 0, Omega_oz, Mx, My, Mz)
        res = computeAcc(0, 0, Omega_oz, OmegaDot_ox, OmegaDot_oy, OmegaDot_oz)
        ad[0][0][t] = res[0][0]
        ad[0][1][t] = res[1][0]
        ad[0][2][t] = res[2][0]

        [OmegaDot_ox, OmegaDot_oy, OmegaDot_oz] = f.getAngularAcc(0, Omega_oy, p0, Mx, My, Mz)
        res = computeAcc(0, Omega_oy, p0, OmegaDot_ox, OmegaDot_oy, OmegaDot_oz)
        ad[1][0][t] = res[0][0]
        ad[1][1][t] = res[1][0]
        ad[1][2][t] = res[2][0]

        [OmegaDot_ox, OmegaDot_oy, OmegaDot_oz] = f.getAngularAcc(Omega_ox, 0, p0, Mx, My, Mz)
        res = computeAcc(Omega_ox, 0, p0, OmegaDot_ox, OmegaDot_oy, OmegaDot_oz)
        ad[2][0][t] = res[0][0]
        ad[2][1][t] = res[1][0]
        ad[2][2][t] = res[2][0]
    
        [OmegaDot_ox, OmegaDot_oy, OmegaDot_oz] = f.getAngularAcc(Omega_ox, Omega_oy, Omega_oz, Mx, My, Mz)
        res = computeAcc(Omega_ox, Omega_oy, Omega_oz, OmegaDot_ox, OmegaDot_oy, OmegaDot_oz)
        ad[3][0][t] = res[0][0]
        ad[3][1][t] = res[1][0]
        ad[3][2][t] = res[2][0]

        OmegaDot_oz = dataOmegaDot[t][1]
        OmegaDot_oy = dataOmegaDot[t][2]
        OmegaDot_ox = dataOmegaDot[t][3]
    
    # Plot results
    f.plotAngVel(dataOmega.T[1], ad[0][:][:], '/Users/sergiocollibars/Desktop/CSML/figures/main3_fig', 'omega_z')
    f.plotAngVel(dataOmega.T[2], ad[1][:][:], '/Users/sergiocollibars/Desktop/CSML/figures/main3_fig', 'omega_y')
    f.plotAngVel(dataOmega.T[3], ad[2][:][:], '/Users/sergiocollibars/Desktop/CSML/figures/main3_fig', 'omega_x')
    f.plotAngVel(np.sqrt(pow(dataOmega.T[3],2) + pow(dataOmega.T[2],2) + pow(dataOmega.T[1],2)), ad[3][:][:], '/Users/sergiocollibars/Desktop/CSML/figures/main3_fig', 'omega_norm')
    
    return 0


# --------------------------------------------------------------------#
#                               MAIN4()                               #
# --------------------------------------------------------------------#

def main4():  

    # Description: Evaluate QGG mesuraments at maximum gamma postion for different types of
    #   torque perturbations. (Constant, gaussian, periodic)

    # orbit parameters
    r_o = np.array([[c.R_E + r], [0],[0]])                      # Orbit radius
    n = math.sqrt(c.G * c.M_E / pow(f.norm(r_o), 3))            # Orbit mean motion

    # Input data
    t_min = 0
    t_max = 6548
    At = 1                                                      # Orbit time step
    deltaT = 1                                                  # Integration time step

    N = int((t_max - t_min) / At)
    time = np.array(np.zeros((1, N)))                           # Time vector

    pert  = 'constant'                                          # perturbation type 
    
    # directories 
    path  = '/Users/sergiocollibars/Desktop/CSML/codes/python_dataFiles/'
    path2 = '/Users/sergiocollibars/Desktop/CSML/figures/main4_fig'

    # For each axis in M
    for k in range(1, 4):
        print(f'Iteration {k}: \n')
        # Define matrix
        M = np.array(np.zeros((4, 1)))  
        omega = np.array(np.zeros((4, 1)))
        omegaDot = np.array(np.zeros((4, 1))) 

        Ma = np.array(np.zeros((4, 1)))  
        Mb = np.array(np.zeros((4, 1)))  

        fa = np.array(np.zeros((4, 1)))  
        fb = np.array(np.zeros((4, 1)))  
            
        # create output files
        f.saveFile(path + 'momentData.txt',    0, 'time,M_x,M_y,M_z')
        f.saveFile(path + 'omegaData.txt',     0, 'time,omega_x,omega_y,omega_z')
        f.saveFile(path + 'omegaDotData.txt',  0, 'time,Omega_x,Omega_y,Omega_z')
        f.saveFile(path + 'adData.txt',        0, 'time,ad_x,ad_y,ad_z')
        f.saveFile(path + 'angleData.txt',     0, 'time,theta,beta,gamma')

        omega1 = np.array(np.zeros((4, 1)))                      # previous omega value
        # for eat time step
        for t in tqdm(range(0, N)):
            # print(f'compute time step: {t}\n')
            if (pert == 'constant'):
                funct = f.torqueConst
            
            # compute perturbations
            M = funct(M, k, t, n)
            M[0][0] = t * At

            Ma = funct(Ma, k, t * At - deltaT, n)
            Ma[0][0] = t * At - deltaT

            Mb = funct(Mb, k, t * At + deltaT, n)
            Mb[0][0] = t * At + deltaT

            # compute angular acceleration
            omegaDot = f.getAngularAcc(omega, omegaDot, M)
            omegaDot[0][0] = t * At

            fa = f.getAngularAcc(omega, fa, Ma)
            fa[0][0] = t * At - deltaT

            fb = f.getAngularAcc(omega, fb, Mb)
            fb[0][0] = t * At + deltaT

            # compute angular velocity
            omega = f.simpsonRule(omega1, omega, omegaDot, fa, fb, deltaT)
            omega[0][0] = t * At

            # update previous angular velocity value
            omega1 = omega

            # Modify spacecraft postion
            gamma = omega[3][0] * t * At
            beta  = omega[2][0] * t * At
            theta = omega[1][0] * t * At
            angle = np.array([[t * At],
                             [theta],
                             [beta],
                             [gamma]])

            # compute sensor measurements
            ad = f.computeAcc(omega, omegaDot, r_o, A1, A2, gamma, beta, theta)

            # save values
            f.saveFile(path + 'momentData.txt', M, '')
            f.saveFile(path + 'omegaData.txt', omega, '')
            f.saveFile(path + 'omegaDotData.txt', omegaDot, '')
            f.saveFile(path + 'adData.txt', ad, '')
            f.saveFile(path + 'angleData.txt', angle, '')

            time[0][t] = At * t


        # read and plot values
        moments = f.readFile(path + 'momentData.txt')
        f.plotVarOrbit(time.T, moments.T, path2, f'Moments_M{k}')

        angVel  = f.readFile(path + 'omegaData.txt')
        f.plotVarOrbit(time.T, angVel.T, path2, f'angularVel_M{k}')

        angAcc  = f.readFile(path + 'omegaDotData.txt')
        f.plotVarOrbit(time.T, angAcc.T, path2, f'angularAcc_M{k}')

        adAcc   = f.readFile(path + 'adData.txt')
        f.plotVarOrbit(time.T, adAcc.T, path2, f'sensorAcc_M{k}')

        angVal  = f.readFile(path + 'angleData.txt')
        f.plotVarOrbit(time.T, angVal.T, path2, f'angle_M{k}')

        print('Next Iteration \n')
        

    return 0


if __name__ == "__main__":
    
    print(f'systema args: {sys.argv}')

    # Print initial configuration
    print(f"\n\nLAUNCH MAIN {mainConfig}\n\n")
    print(f"INITAL PARAMETERS \nNumber of gravity parameters: {gravityPotential}")
    print(f"Longitude: {longitude}")
    print(f"Latitude: {latitude}")
    print(f"Orbit radius (over Earth): {r}")
    print(f"Initial swipe degree: {gamma0}")
    print(f"Final swipe degree: {gamma1}")
    print(f"QGG distance sensor 1: {A1}")
    print(f"QGG distance sensor 2: {A2}\n\n")

    if mainConfig == 1:
        main1()
    elif mainConfig == 2:
        main2()
    elif mainConfig == 3:
        main3()
    elif mainConfig == 4:
        main4()
    else:
        raise Exception('Main configuration not correct. Try other...')
