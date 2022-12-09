#!/usr/local/bin/python
import math
import numpy
from os import listdir
import re
import sys
import lib_process_data as lpd

beta = {16:(6.872), 20:(7.035), 24:(7.192), 36:(7.544), 48:(7.793)}

for file in sys.argv[1:]:
    nt = re.sub(r'^(.*)_', '', file); nt = re.sub(r'(^.*?)\..*', r'\1', nt) ; nt = int(nt)
    nt_half = int(nt/2)
    res = numpy.ndarray((3,nt_half))
    tmp = numpy.loadtxt("./"+file, unpack=True)

    for j in range(nt_half):
        res[0,j] = lpd.tauT_imp[nt][j]
        #res[0,j] = tmp[0,j]
        for i in (1,2):
            res[i,j] =  tmp[i,j] * ((1+0.13771856909427574*6/beta[nt])/(1+0.079*6/beta[nt])) #remove old renorm factor and put it new correct one
            res[i,j] = res[i,j] * lpd.norm_corr(tmp[0,j])  #remove G_norm at the tree level improved distance
            res[i,j] = res[i,j] / lpd.EE_norm_latt_flow[nt][j][2] #improve correlator at the original not improved distance
    numpy.savetxt("EE_2015_Nt"+str(nt)+".dat", numpy.swapaxes(res, 0, 1))
