#import numpy as np
#import scipy.linalg as scl
from numpy import sin,cos,array,transpose,matmul
from numpy.linalg import inv
from scipy.linalg import expm



def flowed_correlator_matrix(action, flow, p1,p2,p3,p4, tau1, tau2, alpha = 1, lam = 1):
# arguments: desired action, definition of flow, four-momentum, flow times, gauge fixing

    # define lattice momenta etc
    def p_hat(p):
        return 2*sin(p/2)
    def c(p):
        return cos(p/2)

    # define Kronecker delta for iteration
    def delta(i,j):
        if i == j:
            return 1
        else:
            return 0

    # turn the four imput momenta into arrays of lattice momenta (and their square)
    ps = array([p_hat(p1),p_hat(p2),p_hat(p3),p_hat(p4)])
    psquared = sum(ps**2)
    cs = array([c(p1),c(p2),c(p3),c(p4)])


    # define kernels of plaquette and rectangles as they will be needed for action kernel
    plaquette_kernel = array([[psquared*delta(i,j)-ps[i]*ps[j] for i in range(4)]
                            for j in range(4)])
    rectangle_kernel = 4*array([[(cs[i]**2*psquared+sum(ps**2*cs**2))*delta(i,j)-(cs[i]**2+cs[j]**2)*ps[i]*ps[j] for i in range(4)]
                            for j in range(4)])
    # additionally, the term for gauge fixing the unflowed propagator is needed
    gf_lambda_kernel = lam*array([[ps[i]*ps[j] for i in range(4)] for j in range(4)])


    # define which action kernel is used
    if action == 'Wilson' or action == 'plaquette': ### if plaquette action
        action_kernel = plaquette_kernel
    elif action == 'rectangle': ### if rectangle action
        action_kernel = (1/8)*rectangle_kernel
    elif action == 'LW' or action == 'Lüscher-Weisz' or action == 'Luscher-Weisz' or action == 'Luescher-Weisz': ### if Lüscher-Weisz action
        action_kernel = (5/3)*plaquette_kernel - (1/12)*rectangle_kernel
    else: ### error message if wrong input for action
        return 'No such action kernel is known. Options are Wilson, plaquette, rectangle, LW, Lüscher-Weisz, Luscher-Weisz or Luescher-Weisz.'

    # add gauge fixing term to action kernel
    action_kernel = action_kernel + gf_lambda_kernel

    # (unflowed) propagator is inverse of (gauge fixed) action kernel
    unfl_prop = inv(action_kernel)



    # define gauge fixing term for flow
    gf_alpha_kernel = alpha*array([[ps[i]*ps[j] for i in range(4)] for j in range(4)])

    # define which kernel to define flow is used
    if flow == 'Wilson' or flow == 'plaquette': ### if Wilson flow
        flow_kernel = plaquette_kernel
    elif flow == 'rectangle': ### if rectangle flow
        flow_kernel = (1/8)*rectangle_kernel
    elif flow == 'LW' or flow == 'Lüscher-Weisz' or flow == 'Luscher-Weisz' or flow == 'Luescher-Weisz': ### if Lüscher-Weisz flow
        flow_kernel = (5/3)*plaquette_kernel - (1/12)*rectangle_kernel
    elif flow == 'Zeuthen': ### if Zeuthen flow
        Zpk = array([[psquared*delta(i,j)*ps[i]**2-ps[i]**3*ps[j] for i in range(4)]
                    for j in range (4)]) ### Zeuthen extra term for plaquette
        Zeuthen_pk = plaquette_kernel - (1/12)*Zpk ### full Zeuthen flow expression for plaquette
        Zrk = 4*array([[(cs[i]**2*psquared+sum(ps**2*cs**2))*ps[i]**2*delta(i,j)-(cs[i]**2+cs[j]**2)*ps[i]**3*ps[j]
                    for i in range (4)]
                        for j in range(4)]) ### Zeuthen extra term for rectangle
        Zeuthen_rk = rectangle_kernel - (1/12)*Zrk ### Zeuthen flowed rectangle kernel
        flow_kernel = (5/3)*Zeuthen_pk - (1/12)*Zeuthen_rk ### full Zeuthen flow kernel
    else: ### error message if wrong input for flow
        return 'Unknown flow kernel. Options are Wilson, plaquette, rectangle, Lüscher-Weisz, LW, Luscher-Weisz, Luescher-Weisz and Zeuthen'

	# add gauge fixing term to flow kernel
    flow_kernel = flow_kernel + gf_alpha_kernel

	# obtain matrix exponential of (gauge fixed) flow kernel for both flow times (=heat kernel)
    exp1 = expm(-tau1*flow_kernel)
    if tau1 == tau2:
        exp2 = exp1
    else:
        exp2 = expm(-tau2*flow_kernel)

    # calculate flowed correlator matrix from dot product of action kernel & heat kernels
    flowed_corr_matrix = matmul(exp1,matmul(unfl_prop,transpose(exp2)))

    # return desired component of flowed propagator
    return flowed_corr_matrix