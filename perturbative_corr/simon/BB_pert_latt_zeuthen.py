#!/usr/local/bin/python3.7 -u

import argparse
import sys

import numpy
import numpy as np
import scipy.linalg

import numba
import concurrent.futures


# first define stuff needed for precise integration (shift, weight)
N_IMP = 4
SH_ONE = .2222726231881051504
SH_TWO = .1799620871697125296
WT_ONE = .6957096902749077147
WT_TWO = 1.3042903097250922853


@numba.njit
def SHIFT(i):
    condition = (i % 4)

    if condition < 2:
        if condition == 0:
            addend = -SH_ONE
        else:
            addend = -SH_TWO
    elif condition == 2:
        addend = SH_TWO
    else:
        addend = SH_ONE

    return i + 0.5 + addend


@numba.njit
def BRANGE(x):
    return x ** 2


@numba.njit
def BJACOBIAN(x):
    return 2 * x


@numba.njit
def WEIGHT(i):
    condition = (i % 4)

    if (condition == 0) or (condition == 3):
        return WT_ONE
    else:
        return WT_TWO


@numba.njit
def flowed_correlator_matrix(action, flow, p1, p2, p3, p4, tau1, tau2, alpha=1, lam=1):
    # arguments: desired action, definition of flow, four-momentum, flow times, gauge fixing

    # define lattice momenta etc
    def p_hat(p):
        return 2 * np.sin(p / 2)

    def c(p):
        return np.cos(p / 2)

    # define Kronecker delta for iteration
    def delta(i, j):
        if i == j:
            return 1
        else:
            return 0

    # turn the four imput momenta into arrays of lattice momenta (and their square)
    ps = np.array([p_hat(p1), p_hat(p2), p_hat(p3), p_hat(p4)])
    psquared = sum(ps ** 2)
    cs = np.array([c(p1), c(p2), c(p3), c(p4)])

    # define kernels of plaquette and rectangles as they will be needed for action kernel
    plaquette_kernel = np.array([[psquared * delta(i, j) - ps[i] * ps[j] for i in range(4)]
                                 for j in range(4)])
    rectangle_kernel = 4 * np.array(
        [[(cs[i] ** 2 * psquared + sum(ps ** 2 * cs ** 2)) * delta(i, j) - (cs[i] ** 2 + cs[j] ** 2) * ps[i] * ps[j] for i in range(4)]
         for j in range(4)])
    # additionally, the term for gauge fixing the unflowed propagator is needed
    gf_lambda_kernel = lam * np.array([[ps[i] * ps[j] for i in range(4)] for j in range(4)])

    # define which action kernel is used
    if action == 'Wilson' or action == 'plaquette':  # if plaquette action
        action_kernel = plaquette_kernel
    elif action == 'rectangle':  # if rectangle action
        action_kernel = (1 / 8) * rectangle_kernel
    elif action == 'LW' or action == 'Luscher-Weisz' or action == 'Luescher-Weisz':  # if Luscher-Weisz action
        action_kernel = (5 / 3) * plaquette_kernel - (1 / 12) * rectangle_kernel
    # add gauge fixing term to action kernel
    action_kernel = action_kernel + gf_lambda_kernel

    # (unflowed) propagator is inverse of (gauge fixed) action kernel
    unfl_prop = np.linalg.inv(action_kernel)

    # define gauge fixing term for flow
    gf_alpha_kernel = alpha * np.array([[ps[i] * ps[j] for i in range(4)] for j in range(4)])

    # define which kernel to define flow is used
    if flow == 'Wilson' or flow == 'plaquette':
        flow_kernel = plaquette_kernel
    elif flow == 'rectangle':
        flow_kernel = (1 / 8) * rectangle_kernel
    elif flow == 'LW' or flow == 'Luscher-Weisz' or flow == 'Luescher-Weisz':  # if Luscher-Weisz flow
        flow_kernel = (5 / 3) * plaquette_kernel - (1 / 12) * rectangle_kernel
    elif flow == 'Zeuthen':
        Zpk = np.array([[psquared * delta(i, j) * ps[i] ** 2 - ps[i] ** 3 * ps[j] for i in range(4)]
                        for j in range(4)])  # Zeuthen extra term for plaquette
        Zeuthen_pk = plaquette_kernel - (1 / 12) * Zpk  # full Zeuthen flow expression for plaquette
        Zrk = 4 * np.array([[(cs[i] ** 2 * psquared + sum(ps ** 2 * cs ** 2)) * ps[i] ** 2 * delta(i, j) - (cs[i] ** 2 + cs[j] ** 2) * ps[i] ** 3 * ps[j]
                             for i in range(4)]
                            for j in range(4)])  # Zeuthen extra term for rectangle
        Zeuthen_rk = rectangle_kernel - (1 / 12) * Zrk  # Zeuthen flowed rectangle kernel
        flow_kernel = (5 / 3) * Zeuthen_pk - (1 / 12) * Zeuthen_rk  # full Zeuthen flow kernel

    # add gauge fixing term to flow kernel
    flow_kernel = flow_kernel + gf_alpha_kernel

    # obtain matrix exponential of (gauge fixed) flow kernel for both flow times (=heat kernel)
    with numba.objmode(answer='float64[:,:]'):  # this is necessary to use unsupported scipy functions inside of a jit-compiled block
        answer = scipy.linalg.expm(-tau1 * flow_kernel)
    exp1 = answer

    if tau1 == tau2:
        exp2 = exp1
    else:
        with numba.objmode(answer='float64[:,:]'):
            answer = scipy.linalg.expm(-tau2 * flow_kernel)
        exp2 = answer

    # calculate flowed correlator matrix from dot product of action kernel & heat kernels
    # use contiguous arrays to suppress numba warnings
    flowed_corr_matrix = np.dot(numpy.ascontiguousarray(exp1), np.dot(numpy.ascontiguousarray(unfl_prop), numpy.ascontiguousarray(np.transpose(exp2))))

    # return desired component of flowed propagator
    return flowed_corr_matrix


class ComputationClass:
    def __init__(self, flow_times, ls_t, N_t, t, ls_space, N_space, factor_space, up, printprogress, nproc):
        self._ls_t = ls_t
        self._N_t = N_t
        self._t = t
        self._ls_space = ls_space
        self._N_space = N_space
        self._factor_space = factor_space
        self._up = up
        self._flow_times = flow_times
        self._nproc = nproc
        self._printprogress = printprogress
        self._result = self.parallelization_wrapper()

    def parallelization_wrapper(self):
        with numba.objmode(answer='float64[:,:]'):
            if self._printprogress:
                print("Started parallel work on these flow times:")
            with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
                result = executor.map(self.pass_argument_wrapper, self._flow_times)
            answer = np.asarray(list(result))
        return answer

    def pass_argument_wrapper(self, tau_f):
        return self.actual_corr_computation(tau_f, self._ls_t, self._N_t, self._t, self._ls_space, self._N_space, self._factor_space, self._up, self._printprogress)

    @staticmethod
    @numba.njit
    def actual_corr_computation(tau_f, ls_t, N_t, t, ls_space, N_space, factor_space, up, printprogress):
        if printprogress:
            with numba.objmode():
                print(r'{0:.6f}'.format(tau_f), end=' ')

        # starting correlator value at specific tau, tau_F that will be increased each step
        value = 0

        # iteration over all time momenta; i is temporal summation index
        for i in ls_t:
            val_t = 0  # value for this specific p_t
            p_t = 2 * np.pi * i / N_t  # temporal momentum
            pt_hat = 2 * np.sin(p_t / 2)  # temporal lattice momentum
            cos_at_t = np.cos(2 * np.pi * i * t)  # factor that will be used later, not needed for spatial integration

            # numerical p_t integration, j is summation index
            for j in ls_space:
                val_x = 0  # value for this specific p_x
                u_x = SHIFT(j) / N_space  # define clever sampling position
                wt_x = WEIGHT(j) * factor_space  # weight these sampling positions cleverly
                p_x = BRANGE(u_x) * up  # shift sampling positions closer to origin (bottom sampling)
                px_hat = 2 * np.sin(p_x / 2)  # x component of lattice momentum
                wt_x *= BJACOBIAN(u_x)  # multiply weight with jacobian to account for bottom sampling

                for k in ls_space:
                    val_y = 0
                    u_y = SHIFT(k) / N_space
                    wt_y = WEIGHT(k) * factor_space
                    p_y = BRANGE(u_y) * up
                    py_hat = 2 * np.sin(p_y / 2)
                    wt_y *= BJACOBIAN(u_y)

                    for l in ls_space:
                        u_z = SHIFT(l) / N_space
                        wt_z = WEIGHT(l) * factor_space
                        p_z = BRANGE(u_z) * up
                        pz_hat = 2 * np.sin(p_z / 2)
                        wt_z *= BJACOBIAN(u_z)

                        # actually calculate the integrand at this p_t,p_x,p_y,p_z
                        
                        # array to iterate over
                        p_vec = np.array([0,px_hat,py_hat,pz_hat])
                        # diagonal term of operator matrix
                        diago_term = np.sum(p_vec**2)*np.diag([0,1,1,1])
                        # non-diagonal term of operator matrix
                        nondiago_term = np.array([[p_vec[mu]*p_vec[nu] for mu in range(4)] for nu in range(4)])
                        # full operator matrix
                        operator_matrix = diago_term - nondiago_term
                        
                        # calculate propagator matrix
                        cm = 16 * flowed_correlator_matrix('Wilson', 'Zeuthen', p_t, p_x, p_y, p_z, tau_f, tau_f)
                        
                        # put operator matrix and propagator together
                        val_at_pos = np.sum(operator_matrix*cm)
                                
                        
                        val_y += val_at_pos * wt_z  # modify value at fixed y at each z

                    val_x += val_y * wt_y  # modify value at fixed x at each y
                val_t += val_x * wt_x  # modify value at fixed t at each y

            # addends i == 0, N_t/2 need to be counted only once, not twice
            if i == 0 or i == int(N_t / 2):
                val_t *= 0.5

            # modify correlator value at each t. Afterwards, it's completely summed & integrated
            value += val_t * cos_at_t
        return value

    def getResult(self):
        return self._result


def get_correlators(temp_seps, flow_times, ls_t, ls_space, N_t, N_space, factor_space, up, printprogress: bool, nproc: int):
    correlators = np.empty((len(flow_times), len(temp_seps)), dtype=np.float64)
    for g, t in enumerate(temp_seps):
        if printprogress:
            print('Temporal separation:', r'{0:.4f}'.format(t))
        tmp = ComputationClass(flow_times, ls_t, N_t, t, ls_space, N_space, factor_space, up, printprogress, nproc)
        results = tmp.getResult()
        for f in range(len(flow_times)):
            correlators[f, g] = results[f]
    return correlators


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--N_t', type=int, required=True, help='Temporal extent of lattice')
    parser.add_argument('--N_space', type=int, default=8, help='1/8 of sampling points that will be used for the momentum integration in each'
                                                               ' direction')
    parser.add_argument('--flowradii_file', type=str, required=True,
                        help='path fo file that contains dimensionless flowradii sqrt(8tauF)T that shall be calculated')
    parser.add_argument('--printprogress', action="store_true")
    parser.add_argument('--outputfile', help='path to outputfile', default='BB_WZnaive.txt', type=str)
    parser.add_argument('--nproc', default=1, type=int, help='max number of concurrently running processes')

    args = parser.parse_args()

    # define flow times
    tmp = np.loadtxt(args.flowradii_file)
    flow_times = args.N_t ** 2 * tmp ** 2 / 8  # transform them into \tau_F/a^2

    # upper limit of spatial integration
    up = np.pi

    args.N_space *= N_IMP  # multiply this number by four

    # list of all indices in temporal addition. Sum runs from 0 to args.N_t/2, will be accounted for by factor 2 if i \neq 0,args.N_t/2
    ls_t = np.asarray(range(int(args.N_t / 2 + 1)))

    # list of all indices in spatial integration
    ls_space = np.asarray(range(args.N_space))

    factor_space = up / args.N_space

    # define temporal separations
    temp_seps = np.asarray([n / args.N_t for n in range(1, int(args.N_t / 2 + 1))])  # in the form \tau T

    correlators = get_correlators(temp_seps, flow_times, ls_t, ls_space, args.N_t, args.N_space, factor_space, up, args.printprogress, args.nproc)

    # overall normalization
    overall_factor = args.N_t ** 3 / (24 * np.pi ** 3)
    correlators *= overall_factor

    # save data in .txt file
    np.savetxt(args.outputfile, correlators, header='naive BB observable, Wilson action, Zeuthen flow')


def save_script_call(add_folder=None):
    """ save full script call in logfile """
    import datetime
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S")
    import os
    file_name = os.path.basename(sys.argv[0])
    folderlist = ("./logs/",)
    if add_folder is not None:
        folderlist = (*folderlist, add_folder)
    for folder in folderlist:
        create_folder(folder)
        with open(folder + "/" + file_name + ".log", 'a') as logfile:
            logfile.write(dt_string + "\n" + " ".join(sys.argv) + '\n')
    return


def print_script_call():
    """ save full script call in logfile """
    import datetime
    now = datetime.datetime.now()
    dt_string = now.strftime("%Y/%m/%d %H:%M:%S")
    print("\n" + dt_string + "\n" + " ".join(sys.argv) + '\n')
    return


def create_folder(*paths):
    from pathlib import Path
    for path in paths:
        Path(path).mkdir(parents=True, exist_ok=True)
    return


if __name__ == '__main__':
    print_script_call()
    main()
    save_script_call()
