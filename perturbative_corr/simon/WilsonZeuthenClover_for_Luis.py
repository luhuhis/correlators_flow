from numpy import sin,cos,array,pi,arange,diag,append,savetxt

import flowed_correlators as fc
from py_integrationheader import *

#########################################################################
###                  parameters to choose                             ###
#########################################################################

# 1/8 of sampling points that will be used for the momentum integration in each direction
N_space = 8 # 8 works fine

# lowest used flow time in units \tau_F/a^2
lowest_tau_F = 0
# largest used flow time in units \tau_F/a^2
largest_tau_F = 1
# flow time step size in units \tau_F/a^2
step_size_tau_F = 0.1

# where to store the .txt file that will be generated
data_path = 'WZCl_action_densities.txt'

# whether to print your progress or not; True or False
printprogress = True


##########################################################################
##########################################################################
##########################################################################



# define some useful values needed later

# upper limit of integration
up = pi

# increase number of integration points by four
N_space *= N_IMP
# list to iterate the integration over
ls_space = range(N_space)
# normalization factor for integration
factor_space = up/N_space

# define list of flow times in units of 
taus = arange(lowest_tau_F, largest_tau_F+0.5*step_size_tau_F, step_size_tau_F)

# overall normalisation
overall_factor = 1/(32*pi**4)

# list where data will be stored
action_densities = array([])

# iteration over all flow times
for t in taus:
    if printprogress:
        print('Flow time: '+str(t)+' ', end = '\n')   # print progress
    value = 0   # starting value that will be increase each step

    # iteration over all time momenta
    for i in ls_space:
        val_t = 0   # value for this specific p_t
        u_t = SHIFT(i)/N_space   # choose smart sampling positions
        wt_t = WEIGHT(i)*factor_space   # weight sampling positions appropriately
        p_t = BRANGE(u_t)*up # temporal momentum, preferred at smaller momenta (bottom sampling)
        pt_dot = sin(p_t) # temporal lattice momentum, needed for observable kernel
        ct_hat = cos(p_t/2) # cosine of momentum, needed for observable kernel
        wt_t *= BJACOBIAN(u_t) # modify weight according to bottom sampling
    
        for j in ls_space:
            val_x = 0
            u_x = SHIFT(j)/N_space
            wt_x = WEIGHT(j)*factor_space
            p_x = BRANGE(u_x)*up
            px_dot = sin(p_x)
            cx_hat = cos(p_x/2)
            wt_x *= BJACOBIAN(u_x)
        
            for k in ls_space:
                val_y = 0
                u_y = SHIFT(k)/N_space
                wt_y = WEIGHT(k)*factor_space
                p_y = BRANGE(u_y)*up
                py_dot = sin(p_y)
                cy_hat = cos(p_y/2)
                wt_y *= BJACOBIAN(u_y)
            
                for l in ls_space:
                    u_z = SHIFT(l)/N_space
                    wt_z = WEIGHT(l)*factor_space
                    p_z = BRANGE(u_z)*up
                    pz_dot = sin(p_z)
                    cz_hat = cos(p_z/2)
                    wt_z *= BJACOBIAN(u_z)
                    
                    # store lattice momenta in lists
                    p_dot_vec = array([pt_dot,px_dot,py_dot,pz_dot])
                    c_hat_vec = array([ct_hat,cx_hat,cy_hat,cz_hat])
                    
                    # observable kernel matrix for these momenta
                    diago_matrix = (pt_dot**2+px_dot**2+py_dot**2+pz_dot**2)*diag(c_hat_vec**2) # part proportional to Kronecker delta
                    transversal_matrix = array([[p_dot_vec[mu]*c_hat_vec[mu]*p_dot_vec[nu]*c_hat_vec[nu] # part not proportional to Kronecker delta
                                                 for mu in range(4)] for nu in range(4)])
                    observable_matrix = diago_matrix - transversal_matrix
                    
                    # propagator matrix for these momenta
                    propagator_matrix = 16*fc.flowed_correlator_matrix('Wilson','Zeuthen', p_t,p_x,p_y,p_z, t,t)
                    
                    # dot observable kernel matrix and propagator matrix together
                    val_at_pos = 0
                    for mu in range(4):
                        for nu in range(4):
                            val_at_pos += observable_matrix[mu,nu]*propagator_matrix[nu,mu]
                    
                    # sum over all p_z
                    val_y += val_at_pos*wt_z 
                
                # sum over all p_y
                val_x += val_y*wt_y
            
            # sum over all p_x
            val_t += val_x*wt_x
        
        # sum over all p_t
        value += val_t*wt_t
    
    # store value in list
    action_densities = append(action_densities,value)

# normalize appropriately
action_densities *= overall_factor

# store results & flow times in .txt file
savetxt(data_path,[taus,action_densities], header = 'Clover observable, Wilson action, Zeuthen flow')