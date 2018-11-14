"""

Particle Filter localization sample

author: Atsushi Sakai (@Atsushi_twi)

"""

# [HRS] Obtained from 
#       https://github.com/AtsushiSakai/PythonRobotics/blob/master/Localization/particle_filter/particle_filter.py
#       and modified.

# [HRS] This demo uses simple Monte Carlo method
#       with markers, which has the problem of 
#       fading weights.

import numpy as np
import math
import matplotlib.pyplot as plt

# Estimation parameter of PF
Q = np.diag([0.1])**2  # range error
#Q = np.diag([0.2])**2  # range error
R = np.diag([1.0, np.deg2rad(40.0)])**2  # input error

# [HRS] May be the standard deviation of errors.
#       square makes square error.
#  Simulation parameter
Qsim = np.diag([0.2])**2
Rsim = np.diag([1.0, np.deg2rad(30.0)])**2

DT = 0.05  # time tick [s]
SIM_TIME = 50.0  # simulation time [s]
MAX_RANGE = 20.0  # maximum observation range

# Particle filter parameter
NP = 200  # Number of Particle
NTh = NP / 2.0  # Number of particle for re-sampling

show_animation = False #True
#show_animation = True


def calc_input():
    v = 1.0  # [m/s]
    yawrate = 0.1  # [rad/s]
    u = np.array([[v, yawrate]]).T
    return u


def observation(xTrue, xd, u, RFID):

    xTrue = motion_model(xTrue, u)

    # add noise to gps x-y
    #z = np.zeros((0, 3))
    z = np.zeros((0, 5))

    for i in range(RFID.shape[0]):

        #dx = xTrue[0, 0] - RFID[i, 0]
        #dy = xTrue[1, 0] - RFID[i, 1]
        #d = math.sqrt(dx**2 + dy**2)
        dv = xTrue[:2,0].T - RFID[i]
        d = np.linalg.norm(dv)
        if d <= MAX_RANGE:
            # [HRS] Why use square of a diagnal matrix here ? What does it mean ?
            dn = d + np.random.randn() * Qsim[0, 0]  # add noise
            zi = np.array([[dn, RFID[i, 0], RFID[i, 1], dv[0], dv[1]]])
            z = np.vstack((z, zi))

    # add noise to input
    ud1 = u[0, 0] + np.random.randn() * Rsim[0, 0]
    ud2 = u[1, 0] + np.random.randn() * Rsim[1, 1]
    ud = np.array([[ud1, ud2]]).T

    xd = motion_model(xd, ud)

    return xTrue, z, xd, ud


def motion_model(x, u):

    # [HRS] What the duck is this ???
    F = np.array([[ 1.0 , 0   , 0   , 0 ],
                  [ 0   , 1.0 , 0   , 0 ],
                  [ 0   , 0   , 1.0 , 0 ],
                  [ 0   , 0   , 0   , 0 ]])

    # [HRS] Why 1.0 at B[3,0] ???
    B = np.array([[ DT * math.cos(x[2, 0]) , 0.0 ],
                  [ DT * math.sin(x[2, 0]) , 0.0 ],
                  [ 0.0                    , DT  ],
                  [ 1.0                    , 0.0 ]])

    x = F.dot(x) + B.dot(u)

    return x


def gauss_likelihood(x, sigma):
    p = 1.0 / math.sqrt(2.0 * math.pi * sigma ** 2) * \
        math.exp(-x ** 2 / (2 * sigma ** 2))

    return p

# [HRS] I understand how it's calculating covariance matrix
#       here now, but why covariance is defined like that?
def calc_covariance(xEst, px, pw):
    cov = np.zeros((3, 3))

    for i in range(px.shape[1]):
        dx = (px[:, i] - xEst)[0:3]
        cov += pw[0, i] * dx.dot(dx.T)

    return cov


def pf_localization(px, pw, xEst, PEst, z, u):
    """
    Localization with Particle filter
    """

    for ip in range(NP):
        x = np.array([px[:, ip]]).T
        w = pw[0, ip]
        #  Predict with random input sampling
        ud1 = u[0, 0] + np.random.randn() * Rsim[0, 0]
        ud2 = u[1, 0] + np.random.randn() * Rsim[1, 1]
        ud = np.array([[ud1, ud2]]).T
        x = motion_model(x, ud)

        #  Calc Inportance Weight
        for i in range(z.shape[0]):
            #dx = x[0, 0] - z[i, 1]
            #dy = x[1, 0] - z[i, 2]
            #prez = math.sqrt(dx**2 + dy**2)
            #dz = prez - z[i, 0]
            dz = np.linalg.norm( x.T[0,:2] - z[i,1:3] ) - z[i,0]
            #dz = np.linalg.norm( ( x.T[0,:2] - z[i,1:3] ) - z[i,3:6] )
            w = w * gauss_likelihood(dz, math.sqrt(Q[0, 0])) # just a guassian bell-shaped function around x=0

        px[:, ip] = x[:, 0]
        pw[0, ip] = w

    pw_sum = pw.sum()
    if(pw_sum):
        pw = pw / pw.sum()  # normalize

    #print pw.shape, pw.sum()
    
    xEst = px.dot(pw.T)
    PEst = calc_covariance(xEst, px, pw)

    px, pw = resampling(px, pw)

    return xEst, PEst, px, pw


def resampling(px, pw):
    """
    low variance re-sampling
    """

    Neff = 1.0 / (pw.dot(pw.T))[0, 0]  # Effective particle number
    if Neff < NTh:
        wcum = np.cumsum(pw)
        base = np.cumsum(pw * 0.0 + 1 / NP) - 1 / NP
        resampleid = base + np.random.rand(base.shape[0]) / NP

        inds = []
        ind = 0
        for ip in range(NP):
            while resampleid[ip] > wcum[ind]:
                ind += 1
            inds.append(ind)

        px = px[:, inds]
        pw = np.zeros((1, NP)) + 1.0 / NP  # init weight

    return px, pw


def plot_covariance_ellipse(xEst, PEst):
    Pxy = PEst[0:2, 0:2]
    #print PEst, Pxy,
    eigval, eigvec = np.linalg.eig(Pxy)
    #print eigval, eigvec

    if eigval[0] >= eigval[1]:
        bigind = 0
        smallind = 1
    else:
        bigind = 1
        smallind = 0

    t = np.arange(0, 2 * math.pi + 0.1, 0.1)

    # eigval[bigind] or eiqval[smallind] were occassionally negative numbers extremely
    # close to 0 (~10^-20), catch these cases and set the respective variable to 0
    try:
        a = math.sqrt(eigval[bigind])
    except ValueError:
        a = 0

    try:
        b = math.sqrt(eigval[smallind])
    except ValueError:
        b = 0

    x = [a * math.cos(it) for it in t]
    y = [b * math.sin(it) for it in t]
    angle = math.atan2(eigvec[bigind, 1], eigvec[bigind, 0])
    R = np.array([[math.cos(angle), math.sin(angle)],
                   [-math.sin(angle), math.cos(angle)]])
    fx = R.dot(np.array([[x, y]]))
    px = np.array(fx[0, :] + xEst[0, 0]).flatten()
    py = np.array(fx[1, :] + xEst[1, 0]).flatten()
    plt.plot(px, py, "--r")


def main():
    print(__file__ + " start!!")

    time = 0.0

    # RFID positions [x, y]
    RFID = np.array([[10.0, 0.0],
                     [10.0, 10.0],
                     [0.0, 15.0],
                     [-5.0, 20.0]])

    # State Vector [x y yaw v]'
    xEst = np.zeros((4, 1))
    xTrue = np.zeros((4, 1))
    PEst = np.eye(4)

    px = np.zeros((4, NP))  # Particle store
    pw = np.zeros((1, NP)) + 1.0 / NP  # Particle weight, mean distribution
    xDR = np.zeros((4, 1))  # Dead reckoning

    # history
    hxEst = xEst
    hxTrue = xTrue
    hxDR = xTrue

    while True:
        #print "1)","="*50
        time += DT
        u = calc_input()

        xTrue, z, xDR, ud = observation(xTrue, xDR, u, RFID)

        if(0):
            print px.shape
            print pw.shape
            print xEst
            print PEst

        xEst, PEst, px, pw = pf_localization(px, pw, xEst, PEst, z, ud)

        #print "2)","-"*50

        if(0):
            print px.shape
            print pw.shape
            print xEst
            print PEst
        

        # store data history
        hxEst = np.hstack((hxEst, xEst))
        hxDR = np.hstack((hxDR, xDR))
        hxTrue = np.hstack((hxTrue, xTrue))

        #print "3)","-"*50

        running = (SIM_TIME >= time)
            
        if (show_animation or 
            (not running)):
            plt.cla()

            for i in range(len(z[:, 0])):
                plt.plot([xTrue[0, 0], z[i, 1]], [xTrue[1, 0], z[i, 2]], "-", color="orange")
            plt.plot(RFID[:, 0], RFID[:, 1], "*", color="orange")
            plt.plot(px[0, :], px[1, :], ".r")
            plt.plot(np.array(hxTrue[0, :]).flatten(),
                     np.array(hxTrue[1, :]).flatten(), "-b")
            plt.plot(np.array(hxDR[0, :]).flatten(),
                     np.array(hxDR[1, :]).flatten(), "-k")
            plt.plot(np.array(hxEst[0, :]).flatten(),
                     np.array(hxEst[1, :]).flatten(), "-r")
            plot_covariance_ellipse(xEst, PEst)
            plt.axis("equal")
            plt.grid(True)
            plt.pause(0.00001)
        if (not running):
            break
    plt.show()

if __name__ == '__main__':
    main()
