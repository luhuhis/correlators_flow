#!/usr/local/bin/python3.7 -u

import argparse
import sys

import numpy as np
import scipy.linalg

import numba
import concurrent.futures


# this is necessary to use unsupported scipy functions inside of a jit-compiled function
@numba.njit(cache=True)
def scipy_linalg_expm_wrapper(arg):
    with numba.objmode(answer='float64[:,:]'):
        answer = scipy.linalg.expm(arg)
    return answer


# this function has to be outside the ComputationClass class, as numba doesn't support calling other functions that contain "with numba.objmode" from
# nested functions
@numba.njit(cache=True)
def flowed_correlator_matrix(gauge_action, flow_action, p1, p2, p3, p4, tau1, tau2, alpha=1, lam=1):
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
    psquared = np.sum(ps ** 2)
    cs = np.array([c(p1), c(p2), c(p3), c(p4)])

    # define kernels of plaquette and rectangles as they will be needed for action kernel
    plaquette_kernel = np.array([[psquared * delta(i, j) - ps[i] * ps[j] for i in range(4)]
                                 for j in range(4)])
    rectangle_kernel = 4 * np.array(
        [[(cs[i] ** 2 * psquared + np.sum(ps ** 2 * cs ** 2)) * delta(i, j) - (cs[i] ** 2 + cs[j] ** 2) * ps[i] * ps[j] for i in range(4)]
         for j in range(4)])
    # additionally, the term for gauge fixing the unflowed propagator is needed
    gf_lambda_kernel = lam * np.array([[ps[i] * ps[j] for i in range(4)] for j in range(4)])

    # define which action kernel is used
    if gauge_action == 'Wilson' or gauge_action == 'plaquette':  # if plaquette action
        action_kernel = plaquette_kernel
    elif gauge_action == 'rectangle':  # if rectangle action
        action_kernel = (1 / 8) * rectangle_kernel
    elif gauge_action == 'LW' or gauge_action == 'Lüscher-Weisz' or gauge_action == 'Luscher-Weisz' or gauge_action == 'Luescher-Weisz':  # if Lüscher-Weisz action
        action_kernel = (5 / 3) * plaquette_kernel - (1 / 12) * rectangle_kernel
    # add gauge fixing term to action kernel
    action_kernel = action_kernel + gf_lambda_kernel

    # (unflowed) propagator is inverse of (gauge fixed) action kernel
    unfl_prop = np.ascontiguousarray(np.linalg.inv(action_kernel))

    # define gauge fixing term for flow
    gf_alpha_kernel = alpha * np.array([[ps[i] * ps[j] for i in range(4)] for j in range(4)])

    # define which kernel to define flow is used
    if flow_action == 'Wilson' or flow_action == 'plaquette':
        flow_kernel = plaquette_kernel
    elif flow_action == 'rectangle':
        flow_kernel = (1 / 8) * rectangle_kernel
    elif flow_action == 'LW' or flow_action == 'Lüscher-Weisz' or flow_action == 'Luscher-Weisz' or flow_action == 'Luescher-Weisz':  # if Lüscher-Weisz flow
        flow_kernel = (5 / 3) * plaquette_kernel - (1 / 12) * rectangle_kernel
    elif flow_action == 'Zeuthen':
        Zpk = np.array([[psquared * delta(i, j) * ps[i] ** 2 - ps[i] ** 3 * ps[j] for i in range(4)]
                        for j in range(4)])  # Zeuthen extra term for plaquette
        Zeuthen_pk = plaquette_kernel - (1 / 12) * Zpk  # full Zeuthen flow expression for plaquette
        Zrk = 4 * np.array(
            [[(cs[i] ** 2 * psquared + np.sum(ps ** 2 * cs ** 2)) * ps[i] ** 2 * delta(i, j) - (cs[i] ** 2 + cs[j] ** 2) * ps[i] ** 3 * ps[j]
              for i in range(4)]
             for j in range(4)])  # Zeuthen extra term for rectangle
        Zeuthen_rk = rectangle_kernel - (1 / 12) * Zrk  # Zeuthen flowed rectangle kernel
        flow_kernel = (5 / 3) * Zeuthen_pk - (1 / 12) * Zeuthen_rk  # full Zeuthen flow kernel

    # add gauge fixing term to flow kernel
    flow_kernel = flow_kernel + gf_alpha_kernel

    # obtain matrix exponential of (gauge fixed) flow kernel for both flow times (=heat kernel)
    exp1 = np.ascontiguousarray(scipy_linalg_expm_wrapper(-tau1 * flow_kernel))

    if tau1 == tau2:
        exp2 = exp1
    else:
        exp2 = np.ascontiguousarray(scipy_linalg_expm_wrapper(-tau2 * flow_kernel))

    # calculate flowed correlator matrix from dot product of action kernel & heat kernels
    flowed_corr_matrix = np.dot(exp1, np.dot(unfl_prop, np.transpose(exp2)))

    # return desired component of flowed propagator
    return flowed_corr_matrix


class ComputationClass:
    def __init__(self, flow_times, gauge_action, flow_action, corr, N_t, N_space, printprogress, nproc):
        self._gauge_action = gauge_action
        self._flow_action = flow_action
        self._corr = corr
        self._N_t = N_t
        self._N_space = N_space
        self._flow_times = flow_times
        self._nproc = nproc
        self._printprogress = printprogress
        self._result = self.parallelization_wrapper()

    def parallelization_wrapper(self):
        if self._printprogress:
            print("Started parallel work on these flow times (tau_F/a^2):")
        with concurrent.futures.ProcessPoolExecutor(max_workers=self._nproc) as executor:
            result = executor.map(self.pass_argument_wrapper, self._flow_times)
        print("")
        answer = np.asarray(list(result))
        return answer

    def pass_argument_wrapper(self, tau_f):
        if self._printprogress:
            print(r'{0:.6f}'.format(tau_f), end=' ')
        return self.actual_corr_computation(tau_f, self._gauge_action, self._flow_action, self._corr, self._N_t, self._N_space)

    @staticmethod
    @numba.njit(cache=True)
    def actual_corr_computation(tau_f, gauge_action, flow_action, corr, N_t, N_space):

        # first define stuff needed for precise integration (shift, weight)
        SH_ONE = .2222726231881051504
        SH_TWO = .1799620871697125296
        WT_ONE = .6957096902749077147
        WT_TWO = 1.3042903097250922853

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

        def BRANGE(x):
            return x ** 2

        def BJACOBIAN(x):
            return 2 * x

        def WEIGHT(i):
            condition = (i % 4)
            if (condition == 0) or (condition == 3):
                return WT_ONE
            else:
                return WT_TWO

        N_space *= 4  # the integration quadrature requires 4 intervals, so we need to make sure this is divisble by 4
        up = np.pi  # upper limit of spatial integration
        factor_space = up / N_space

        # starting correlator value at specific tau, tau_F that will be increased each step
        value = np.zeros(N_t)
        temp_seps = np.asarray([n / N_t for n in range(1, int(N_t / 2 + 1))])  # in the form \tau T

        # list of all indices in temporal addition. Sum runs from 0 to args.N_t/2, will be accounted for by factor 2 if i \neq 0,args.N_t/2
        ls_t = range(int(N_t / 2 + 1))
        ls_space = range(N_space)  # list of all indices in spatial integration

        # iteration over all time momenta; i is temporal summation index
        for i in ls_t:
            val_t = 0  # value for this specific p_t
            p_t = 2 * np.pi * i / N_t  # temporal momentum
            pt_hat = 2 * np.sin(p_t / 2)  # temporal lattice momentum

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

                        if corr == "EE":
                            # actually calculate the integrand at this p_t,p_x,p_y,p_z
                            # calculate correlator matrix
                            cm = 16 * flowed_correlator_matrix(gauge_action, flow_action, p_t, p_x, p_y, p_z, tau_f, tau_f)
                            # first term of EE correlator
                            term1 = (px_hat ** 2 + py_hat ** 2 + pz_hat ** 2) * cm[0, 0] + pt_hat ** 2 * (cm[1, 1] + cm[2, 2] + cm[3, 3])
                            # second term of EE correlator
                            term2 = pt_hat * (px_hat * (cm[0, 1] + cm[1, 0]) + py_hat * (cm[0, 2] + cm[2, 0]) + pz_hat * (cm[0, 3] + cm[3, 0]))
                            val_at_pos = term1 - term2  # difference of first and second term
                        elif corr == "BB":
                            # actually calculate the integrand at this p_t,p_x,p_y,p_z

                            # array to iterate over
                            p_vec = np.array([0, px_hat, py_hat, pz_hat])

                            # diagonal term of operator matrix
                            pvecsum = np.sum(p_vec ** 2)
                            diago_term = np.diag(np.array([0, pvecsum, pvecsum, pvecsum]))

                            # non-diagonal term of operator matrix
                            nondiago_term = np.zeros((4, 4))
                            for nu in range(4):
                                for mu in range(4):
                                    nondiago_term[mu, nu] = p_vec[mu] * p_vec[nu]

                            # full operator matrix
                            operator_matrix = diago_term - nondiago_term

                            # calculate propagator matrix
                            cm = 16 * flowed_correlator_matrix(gauge_action, flow_action, p_t, p_x, p_y, p_z, tau_f, tau_f)

                            # put operator matrix and propagator together
                            val_at_pos = np.sum(operator_matrix * cm)

                        val_y += val_at_pos * wt_z  # modify value at fixed y at each z

                    val_x += val_y * wt_y  # modify value at fixed x at each y
                val_t += val_x * wt_x  # modify value at fixed t at each y

            # addends i == 0, N_t/2 need to be counted only once, not twice
            if i == 0 or i == int(N_t / 2):
                val_t *= 0.5

            # modify correlator value at each tauT
            for m, t in enumerate(temp_seps):
                cos_at_t = np.cos(2 * np.pi * i * t)
                value[m] += val_t * cos_at_t

        return value

    def getResult(self):
        return np.asarray(self._result)


def get_correlators(flowtimes_file: str, gauge_action: str, flow_action: str, corr: str, N_t: int, N_space: int, printprogress: bool, nproc: int):
    flow_times = np.loadtxt(flowtimes_file)
    correlators = np.empty((len(flow_times), int(N_t/2)), dtype=np.float64)
    tmp = ComputationClass(flow_times, gauge_action, flow_action, corr, N_t, N_space, printprogress, nproc)
    results = tmp.getResult()
    for g in range(int(N_t/2)):
        for f in range(len(flow_times)):
            correlators[f, g] = results[f, g]

    if corr == "EE":
        sign = -1
    elif corr == "BB":
        sign = 1

    # overall normalization
    overall_factor = sign * N_t ** 3 / (24 * np.pi ** 3)
    correlators *= overall_factor

    return correlators


def add_args(parser):
    """ arguments that are relevant for get_correlators() """
    parser.add_argument('--Nt', type=int, required=True, help='Temporal extent of lattice')
    parser.add_argument('--Nspace', type=int, default=8, help='1/8 of sampling points that will be used for the momentum integration in each'
                                                               ' direction. default is 8 (using 16 changes results only after the ~6th decimal place).')
    parser.add_argument('--flowtimes_file', type=str, required=True,
                        help='path fo file that contains dimensionless flowtimes tauF/a^2 that shall be calculated')
    parser.add_argument('--printprogress', action="store_true")

    parser.add_argument('--nproc', default=1, type=int, help='max number of concurrently running processes')
    parser.add_argument('--flow_action', choices=['Wilson', 'Zeuthen', 'rectangle', 'LW'], type=str, required=True)
    parser.add_argument('--gauge_action', choices=['Wilson', 'rectangle', 'LW'], type=str, required=True)
    parser.add_argument('--corr', choices=['EE', 'BB'], type=str, required=True)


def main():

    # parse arguments
    parser = argparse.ArgumentParser()
    add_args(parser)
    parser.add_argument('--outputpath', help='path to output folder', default='./', type=str, required=True)
    args = parser.parse_args()

    # compute correlators
    correlators = get_correlators(args.flowtimes_file, args.gauge_action, args.flow_action, args.corr, args.Nt, args.Nspace, args.printprogress, args.nproc)

    # save data in text file
    create_folder(args.outputpath)
    np.savetxt(args.outputpath+"/"+args.corr+"_pert_latt_"+args.flow_action+"_"+args.gauge_action+"_Nt"+str(args.Nt)+".dat", correlators, header='rows are flow times, columns are temp seps')


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
    """ print full script call to std out """
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
