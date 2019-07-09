#!/usr/bin/python3
import math
import numpy

def free_corr( tauT ):
    norm = math.pi**2 * ( math.cos(math.pi*tauT)**2 / math.sin(math.pi*tauT)**4 + 1/(3*math.sin(math.pi*tauT)**2))
    return norm

def normalize_data( data, beta ):
    datashape = data[0].shape
    Ntau = datashape[1]*2
    for k in range(0,len(data)):
        for i in range(0, datashape[0]):
            for j in range(0, datashape[1]):
                data[k][i,j] = data[k][i,j] / free_corr((j+1)/Ntau) * Ntau**4 * (1+0.079*6/beta) 
                #if i == 0 : data[i,j] *= (1+0.079*6/beta) 

def normalize_EE( data, beta ):
    datashape = data.shape
    Ntau = datashape[1]*2
    for i in range(0, datashape[0]):
        for j in range(0, datashape[1]):
            data[i,j] = data[i,j] / free_corr((j+1)/Ntau) * Ntau**4 * (1+0.079*6/beta) 
            #if i == 0 : data[i,j] *= (1+0.079*6/beta) 
def get_mean_and_jackknife( data, noBlocks ):
    ndata=len(data)
    blocksize = ndata / noBlocks
    if not blocksize.is_integer():
        exit("Error in jackknife: len(data)/noBlocks is not an integer")
    blocksize = int(blocksize)
    datashape = data[0].shape
    subset = {}

    """calc mean"""
    mean = numpy.zeros(datashape)
    for i in range(0,ndata):
        mean += data[i]
    mean /= ndata

    """calc subset"""
    for i in range(0,noBlocks):
        subset[i] = numpy.zeros(datashape)
        for j in range(0,ndata-blocksize):
            subset[i] += data[(i*blocksize+j)%ndata]
        subset[i] /= ndata-blocksize

    """calc error"""
    err = numpy.zeros(datashape)
    for i in range(0,noBlocks):
        err += (mean-subset[i])*(mean-subset[i])
    err = numpy.sqrt(err/noBlocks*(noBlocks-1))

    return mean, err

def jackknife_EE( numerator, polyakov, noBlocks ):
    n_numerator=len(numerator)
    n_polyakov=len(polyakov)
    if n_numerator != n_polyakov :
        exit("Error in jackknife: n_numerator != n_polyakov")
    
    ndata = n_numerator
    blocksize = ndata / noBlocks
    if not blocksize.is_integer():
        exit("Error in jackknife: len(data)/noBlocks is not an integer")
    blocksize = int(blocksize)
    
    numerator_shape = numerator[0].shape
    polyakov_shape = polyakov[0].shape
    
    numerator_subset = {}
    polyakov_subset = {}
    corr_subset = {}

    """calc mean"""
    numerator_mean = numpy.zeros(numerator_shape)
    polyakov_mean = numpy.zeros(polyakov_shape)
    for i in range(0,ndata):
        numerator_mean += numerator[i]
        polyakov_mean += polyakov[i]
    for i in range(0,polyakov_shape[0]):
        numerator_mean[i,:] /= polyakov_mean[i]
    corr_mean = numerator_mean

    """calc subset"""
    for i in range(0,noBlocks):
        numerator_subset[i] = numpy.zeros(numerator_shape)
        polyakov_subset[i] = numpy.zeros(polyakov_shape)
        for j in range(0,ndata-blocksize):
            numerator_subset[i] += numerator[(i*blocksize+j)%ndata]
            polyakov_subset[i] += polyakov[(i*blocksize+j)%ndata]
        for k in range(0,polyakov_shape[0]):
            numerator_subset[i][k,:] /= polyakov_subset[i][k]
        corr_subset[i] = numerator_subset[i]

    """calc error"""
    corr_err = numpy.zeros(numerator_shape)
    for i in range(0,noBlocks):
        corr_err += (corr_mean-corr_subset[i])**2
    corr_err = numpy.sqrt(corr_err/noBlocks*(noBlocks-1))

    return corr_mean, corr_err
