#---------------------------------------------------------------------#
#                          FUNCTIONS                                  # 
#                                                                     #
# Author: Sergio Coll                                                 #
# Date: 08/25/22                                                      #
# Description: Set of functions used in main.py                       #
#---------------------------------------------------------------------#
from cmath import inf
import constants as c
import math 
import numpy as np
import matplotlib.pyplot as plt

# Gets potential gravity for a given orbit elements.
# Uses Legendre funtion.
def getPotential(terms_num, G, M, r, Lambda, Phi):
    # r: Orbit radius
    # terms_num: Number of terms in Legendre functions
    # G: Gravity constant
    # M: Planet mass
    # Lamda: Orbit longitude from inertial frame from planet centre
    # Phi: Orbit latitude from inertial frame from planet centre
    
    U = G * M / r

    if (terms_num > 1):
        raise Exception(" Parameters for Legendre function not implemented yet ...")

    return U

# Computes angular velocity for a circular orbit around a planet
def getAngularVel(G, M, r):
    # r: Orbit radius, module
    # terms_num: Number of terms in Legendre functions
    # G: Gravity constant
    # M: Planet mass
    
    Omega_x = 0
    Omega_y = 0
    Omega_z = math.sqrt(G * M / pow(r, 3))

    return [Omega_x, Omega_y, Omega_z]

# Computes angular acceleration using linearized EOM in a satellite
def getAngularAcc(omega, omegaDot, M):
    # omega : Angular velocity vector [rad / s]
    # M: Torque vector [Nm]
    # omegaDot: Angular acceleration vector [rad / s^2]

    Ix = 1000                              # Satellite initertia x axis [kg m^3]
    Iy = 1000                             # Satellite initertia y axis [kg m^3]
    Iz = 1000                              # Satellite initertia z axis [kg m^3]
    
    # omegaDot[1][0] = omega[2][0] * (Iy - Iz) * omega[3][0] / Ix + M[1][0] / Ix
    # omegaDot[2][0] = omega[1][0] * (Iz - Ix) * omega[3][0] / Iy + M[2][0] / Iy
    # omegaDot[3][0] = omega[1][0] * (Ix - Iy) * omega[2][0] / Iz + M[3][0] / Iz

    omegaDot[1][0] =  M[1][0] / Ix
    omegaDot[2][0] =  M[2][0] / Iy
    omegaDot[3][0] =  M[3][0] / Iz

    return omegaDot

# Computes the sensor measurement
def computeAcc(omega, Omega, r, A1, A2, gamma, beta, theta):
    # omega: Angular velocity vector [rad / s]
    # Omega: Angular acceleration vector [rad / s^2]
    # r: Orbit radius   [m]
    # A1: Distance vector form COM for sensor 1 [m]
    # A2: Distance vector from COM for sensor 2 [m]
    # gamma: Angle from Xs axis [rad]

    # sensor measurment
    ad = np.array(np.zeros((4, 1)))

    # Orbit angular velocity
    p0 = math.sqrt(c.G * c.M_E / pow(norm(r), 3))

    # Compute gravity gradient tensor
    U = np.array([[2 * c.G * c.M_E / pow(norm(r), 3), 0, 0],
            [0, -c.G * c.M_E / pow(norm(r), 3), 0],
            [0, 0, -c.G * c.M_E / pow(norm(r), 3)]])

    # Construct angular acceleration matrixs
    Omega2 = np.array([[-pow(omega[3][0] + p0,2) -pow(omega[2][0],2), omega[1][0] * omega[2][0], omega[1][0] * (omega[3][0] + p0)],
                [omega[1][0] * omega[2][0], -pow(omega[3][0] + p0,2) -pow(omega[1][0],2), omega[2][0] * (omega[3][0] + p0)],
                [omega[1][0] * (omega[3][0] + p0), omega[2][0] * (omega[3][0] + p0), -pow(omega[1][0],2) -pow(omega[2][0],2)]])
    
    OmegaDot = np.array([[0, -Omega[3][0], Omega[2][0]],
                    [Omega[3][0], 0, -Omega[1][0]],
                    [-Omega[2][0], Omega[1][0], 0]])

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

    R_inv = np.linalg.inv(np.dot(R_x, np.dot(R_y, R_z)))
    
    Delta_r1 = np.dot(R_inv, A1.T)                                            # Sensor 2 radius, Axis: {xs, ys, zs}
    Delta_r2 = np.dot(R_inv, A2.T)                                           # Sensor 2 radius, Axis: {xs, ys, zs}

    # Compute acceleration
    a1 = np.dot(-(U - Omega2 - OmegaDot), Delta_r1)
    a2 = np.dot(-(U - Omega2 - OmegaDot), Delta_r2) 

    ad[0][0] = omega[0][0]                                                    # time [s]
    ad[1][0] = (a1 - a2)[0][0]                                                # ax [m / s^2]
    ad[2][0] = (a1 - a2)[1][0]                                                # ay [m / s^2]
    ad[3][0] = (a1 - a2)[2][0]                                                # ax [m / s^2]

    return ad # x, y, z

# Computes partial derivative of a vector respect to an axis
def partial(f1, f2, deltaPos):
    # f1: vector at point 1
    # f2: vector at point 2
    # deltaPos: Increment in position.
    # Euler: df/dx = (f2 - f1) / (deltaPos)
    
    return (f2 - f1) / deltaPos

# Check infinity values inside matric
def checkInf(f):
    # f: matrix
    
    for i in range(0, len(f)):
        if(f[i] == inf or f[i] == -inf): f[i] = 0
    return f

# Computes norm of a vector
def norm(r):
    # r: vector n components. Column vector [n x 1]
    
    val = 0
    for i in range(len(r)):
        val = pow(r[i][0], 2) + val

    norm = np.sqrt(val)

    return norm

# Integrates angular acceleration using Simpson's rule
def simpsonRule(omega1, omega, omegaDot, fa, fb, deltaT):
    # omega1: previous value of angular velocity [rad / s]
    # omega: angular velocity [rad / s]
    # omegaDot: angular acceleration [rad / s^2]
    # fa: function at point a: t - deltaT
    # fb: function at point b: t + deltaT
    # deltaT: time increment

    omega[1][0] = 2 * deltaT / 6 * (fa[1][0] + 4 * omegaDot[1][0] + fb[1][0]) + omega1[1][0]
    omega[2][0] = 2 * deltaT / 6 * (fa[2][0] + 4 * omegaDot[2][0] + fb[2][0]) + omega1[2][0]
    omega[3][0] = 2 * deltaT / 6 * (fa[3][0] + 4 * omegaDot[3][0] + fb[3][0]) + omega1[3][0]

    return omega

# Computes constant torque
def torqueConst(M, k, t, n):
    # M: torque vector [t, x, y, z]
    # k: axis index [1, 2, 3]
    # t: current time step
    # n: mean orbit motion

    pert = 1000 * pow(n, 2) / (math.pi * 720)         # Perturbation [Nm]
    M[k][0] = pert                                    # Constant torque

    return M

# Plot angular acceleration
def plotAcc(angle, ad, location):
    fig1, axs = plt.subplots(3, 1)
    for j in range(3):

        axs[j].plot(np.rad2deg(angle), ad[j][:])
        axs[j].set_xlim(min(np.rad2deg(angle)), max(np.rad2deg(angle)))
        axs[j].set_xlabel('Angle $\gamma$[deg]')
        axs[j].grid(True)

    axs[0].set_ylabel('$a_{dx}$ [m/s2]')
    axs[1].set_ylabel('$a_{dy}$ [m/s2]')
    axs[2].set_ylabel('$a_{dz}$ [m/s2]')

    plt.savefig(f'{location}/acceleration_components.png')

    fig2, ax = plt.subplots()
    ax.plot(np.rad2deg(angle), np.sqrt(pow(ad[0][:],2) + pow(ad[1][:], 2) + pow(ad[2][:], 2)))
    ax.set_ylabel('Acceleration norm [m/s2]')
    ax.set_xlabel('Angle $\gamma$ [deg]')
    ax.set_title('QGG acceleration |$a_d$|')
    ax.grid(True)

    plt.savefig(f'{location}/acceleration_module.png')

    plt.show()

# Plot angular velocity
def plotAngVel(xAxis, ad, location, name):


    fig2, ax = plt.subplots()
    ax.plot(xAxis[:], np.sqrt(pow(ad[0][:],2) + pow(ad[1][:], 2) + pow(ad[2][:], 2)))
    ax.set_ylabel('Acceleration norm [m/s2]')
    ax.set_xlabel(f'{name}')
    ax.set_title('QGG acceleration |$a_d$|')
    ax.grid(True)

    plt.savefig(f'{location}/{name}.jpg')

    plt.show()

# plot variables along orbit
def plotVarOrbit(xAxis, var, location, name):
    # xAxis: x axis variable
    # var: y axis variable
    # location: path to save file
    # name: plot and file name

    fig2, ax = plt.subplots(4, 1)
    ax[0].plot(xAxis[:], np.sqrt(pow(var[1][:],2) + pow(var[2][:], 2) + pow(var[3][:], 2)))
    ax[0].set_ylabel('Variable norm')
    ax[0].set_xlabel('Time [s]')
    ax[0].set_title(f'{name}')
    ax[0].grid(True)

    ax[1].plot(xAxis[:], var[1][:])
    ax[1].set_ylabel('varaible x axis')
    ax[1].set_xlabel('Time [s]')
    ax[1].grid(True)

    ax[2].plot(xAxis[:], var[2][:])
    ax[2].set_ylabel('varaible y axis')
    ax[2].set_xlabel('Time [s]')
    ax[2].grid(True)

    ax[3].plot(xAxis[:], var[3][:])
    ax[3].set_ylabel('varaible z axis')
    ax[3].set_xlabel('Time [s]')
    ax[3].grid(True)

    fig2.set_size_inches(14.09, 9.68)
    plt.savefig(f'{location}/{name}.jpg')

    # plt.show()
    return 0

# Read txt file
def readFile(path):
    # read txt file. Data is column stack and separeted by ','
    file1 = open(path, 'r')
    Lines = file1.readlines()

    # get number of columns
    colNum = len(Lines[1].split(','))
    data = np.array(np.zeros((len(Lines)-1, colNum)))

    # Fill data matrix
    cont = 0
    for line in Lines:
        a = line.split(',')

        if (a[0] == 'time'): continue
        for i in range(colNum):
            data[cont][i] = float(a[i].replace('\n',''))
        cont = cont + 1
    return data

# function to save data
def saveFile(path, var, header):
    # open txt file. Data is column stack and separeted by ','
    if(header != ''):
        file1 = open(path, 'w')
        file1.write(header + '\n')
        file1.close
    else: 
        file1 = open(path, 'a')
        line = str(format(var[0][0], "10.7E")) + ',' + str(format(var[1][0], "10.7E")) + ',' + str(format(var[2][0], "10.7E")) + ',' + str(format(var[3][0], "10.7E")) + '\n'
        file1.write(line)
        file1.close

    return 0