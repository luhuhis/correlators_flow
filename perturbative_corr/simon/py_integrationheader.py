# first define stuff needed for precise integration (shift, weight)
N_IMP = 4
SH_ONE = .2222726231881051504
SH_TWO = .1799620871697125296
WT_ONE = .6957096902749077147
WT_TWO = 1.3042903097250922853

def SHIFT(i):

    condition = (i%4)

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


def WEIGHT(i):

    condition = (i%4)

    if (condition == 0) or (condition == 3):
        return WT_ONE
    else:
        return WT_TWO



# now functions to rescale the integration range
def SRANGE(x):
    return x**2*(3-2*x)
def SJACOBIAN(x):
    return 6*x*(1-x)

def TRANGE(x):
    return x*(2-x)
def TJACOBIAN(x):
    return 2*(1-x)

def BRANGE(x):
    return x**2
def BJACOBIAN(x):
    return 2*x


# combine the above to define a magnificient integrator
def integrator(function, upper_limit, N, scaling = 'uniform'):

    N *= N_IMP
    ls = range(N)
    value = 0
    i = 0

    if scaling == 'uniform':
        factor = upper_limit/N
        for i in ls:
            x = SHIFT(i)*factor
            wt = WEIGHT(i)*factor
            val_at_x = function(x)
            value += val_at_x*wt
    elif scaling == 'bottom':
        for i in ls:
            y = SHIFT(i)/N
            wt = WEIGHT(i)*upper_limit/N
            x = BRANGE(y)*upper_limit
            wt *= BJACOBIAN(y)
            val_at_x = function(x)
            value += val_at_x*wt
    elif scaling == 'both':
        for i in ls:
            y = SHIFT(i)/N
            wt = WEIGHT(i)*upper_limit/N
            x = SRANGE(y)*upper_limit
            wt *= SJACOBIAN(y)
            val_at_x = function(x)
            value += val_at_x*wt
    elif scaling == 'top':
        for i in ls:
            y = SHIFT(i)/N
            wt = WEIGHT(i)*upper_limit/N
            x = TRANGE(y)*upper_limit
            wt *= TJACOBIAN(y)
            val_at_x = function(x)
            value += val_at_x*wt
    else:
        return 'There is no rescaling of this name. Try uniform, bottom, top, or both.'

    return value
