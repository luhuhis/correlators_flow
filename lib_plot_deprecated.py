import numpy

#TODO: REMOVE THIS SHIT
if False:
    EE_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_s064t16_b0687361.dat")
    EE_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_s080t20_b0703500.dat")
    EE_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_s096t24_b0719200.dat")
    EE_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_s120t30_b0739400.dat")
    EE_36 = numpy.loadtxt(inputfolder+"/s144t36_b0754400/EE_s144t36_b0754400.dat")
    EE = [EE_16, EE_20, EE_24, EE_30, EE_36]
    EE_err_16 = numpy.loadtxt(inputfolder+"/s064t16_b0687361/EE_err_s064t16_b0687361.dat")
    EE_err_20 = numpy.loadtxt(inputfolder+"/s080t20_b0703500/EE_err_s080t20_b0703500.dat")
    EE_err_24 = numpy.loadtxt(inputfolder+"/s096t24_b0719200/EE_err_s096t24_b0719200.dat")
    EE_err_30 = numpy.loadtxt(inputfolder+"/s120t30_b0739400/EE_err_s120t30_b0739400.dat")
    EE_err_36 = numpy.loadtxt(inputfolder+"/s144t36_b0754400/EE_err_s144t36_b0754400.dat")
    EE_err = [EE_err_16, EE_err_20, EE_err_24, EE_err_30, EE_err_36]
    #tauT_30_ext = (*tauT_30,0.5)
    tauT_30_ext = (0.229167, 0.250000, 0.270833, 0.291667, 0.312500, 0.333333, 0.354167, 0.375000, 0.395833, 0.416667, 0.437500, 0.458333, 0.479167, 0.500000)
    flow_radius = numpy.loadtxt(inputfolder+"/s064t16_b0687361/flowradii_s064t16_b0687361.dat")
