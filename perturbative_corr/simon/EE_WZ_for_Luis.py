from numpy import sin,cos,array,pi,append,savetxt,arange,delete

import flowed_correlators as fc
from py_integrationheader import *

#########################################################################
###                  parameters to choose                             ###
#########################################################################

# size of lattice in time direction
N_t = 8

# 1/8 of sampling points that will be used for the momentum integration in each direction
N_space = 8 # 8 works fine

# lowest used flow time in units \sqrt{8\tau_F}T
lowest_tau_F = 0
# largest used flow time in units \sqrt{8\tau_F}T
largest_tau_F = 0.01
# flow time step size in units \sqrt{8\tau_F}T
step_size_tau_F = 0.01

# where to store the .txt file that will be generated
data_path = 'EE_WZnaive.txt'

# whether to print your progress or not; True or False
printprogress = True


##########################################################################
##########################################################################
##########################################################################




# upper limit of spatial integration
up = pi

N_space *= N_IMP # multiply this number by four

# list of all indices in temporal addition. Sum runs from 0 to N_t/2, will be accounted for by factor 2 if i \neq 0,N_t/2
ls_t = range(int(N_t/2+1))

# list of all indices in spatial integration
ls_space = range(N_space)

factor_space = up/N_space

# define temporal separations
temp_seps = array([n/N_t for n in range(1, int(N_t/2+1))]) # in the form \tau T

# define flow times
flow_times = arange(lowest_tau_F,largest_tau_F+0.5*step_size_tau_F,step_size_tau_F) # in the form \sqrt(8\tau_F)T
flow_times = N_t**2*flow_times**2/8 # transform them into \tau_F/a^2

# overall normalization
overall_factor = -N_t**3/(24*pi**3)

# create 2D array that will store correlators at all temporal separations and flow times
correlators = array([0]*len(flow_times)) # force array to be 2D, zero line will be deleted later
correlators = array([correlators])       #          "

# iteration over temporal separations
for t in temp_seps:
    # create list that stores values for different tau_F at fixed tau
    correlator_at_fixed_sep = array([])

    # iteration over all flow times
    for tau_f in flow_times:
        # print progress
        if printprogress:
            print(f'Temporal separation: {t}, Flow time: {tau_f}')
        # starting correlator value at specific tau, tau_F that will be increased each step
        value = 0

        # iteration over all time momenta; i is temporal summation index
        for i in ls_t:
            val_t = 0 # value for this specific p_t
            p_t = 2*pi*i/N_t # temporal momentum
            pt_hat = 2*sin(p_t/2) # temporal lattice momentum
            cos_at_t = cos(2*pi*i*t) # factor that will be used later, not needed for spatial integration

            # numerical p_t integration, j is summation index
            for j in ls_space:
                val_x = 0 # value for this specific p_x
                u_x = SHIFT(j)/N_space # define clever sampling position
                wt_x = WEIGHT(j)*factor_space # weight these sampling positions cleverly
                p_x = BRANGE(u_x)*up # shift sampling positions closer to origin (bottom sampling)
                px_hat = 2*sin(p_x/2) # x component of lattice momentum
                wt_x *= BJACOBIAN(u_x) # multiply weight with jacobian to account for bottom sampling

                for k in ls_space:
                    val_y = 0
                    u_y = SHIFT(k)/N_space
                    wt_y = WEIGHT(k)*factor_space
                    p_y = BRANGE(u_y)*up
                    py_hat = 2*sin(p_y/2)
                    wt_y *= BJACOBIAN(u_y)

                    for l in ls_space:
                        u_z = SHIFT(l)/N_space
                        wt_z = WEIGHT(l)*factor_space
                        p_z = BRANGE(u_z)*up
                        pz_hat = 2*sin(p_z/2)
                        wt_z *= BJACOBIAN(u_z)

                        # actually calculate the integrand at this p_t,p_x,p_y,p_z
                        # calculate correlator matrix
                        cm = 16*fc.flowed_correlator_matrix('Wilson','Zeuthen', p_t,p_x,p_y,p_z, tau_f,tau_f)
                        # first term of EE correlator
                        term1 = (px_hat**2+py_hat**2+pz_hat**2)*cm[0,0] + pt_hat**2*(cm[1,1]+cm[2,2]+cm[3,3])
                        # second term of EE correlator
                        term2 = pt_hat*(px_hat*(cm[0,1]+cm[1,0])+py_hat*(cm[0,2]+cm[2,0])+pz_hat*(cm[0,3]+cm[3,0]))
                        val_at_pos = term1-term2 # difference of first and second term

                        val_y += val_at_pos*wt_z # modify value at fixed y at each z

                    val_x += val_y*wt_y # modify value at fixed x at each y

                val_t += val_x*wt_x # modify value at fixed t at each y

            # addends i == 0, N_t/2 need to be counted only once, not twice
            if i == 0 or i == int(N_t/2):
                val_t *= 0.5

            # modify correlator value at each t. Afterwards, it's completely summed & integrated
            value += val_t*cos_at_t

        # store list at fixed separation with different flow times
        correlator_at_fixed_sep = append(correlator_at_fixed_sep,value)

    # store 2D array containing correlators for all separations and flow times
    correlator_at_fixed_sep = array([correlator_at_fixed_sep])
    correlators = append(correlators, correlator_at_fixed_sep, axis = 0)

# delete the zero line at beginning of 'correlators', normalize appropriately & save G_LO/T^4
correlators = delete(correlators, 0,0)*overall_factor

# save data in .txt file
savetxt(data_path,correlators,
        header = 'naive EE observable, Wilson action, Zeuthen flow')