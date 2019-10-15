#!/usr/bin/python3
import math
import numpy
from os import listdir 

chi_dofs = numpy.loadtxt("/home/altenkort/master/work/data_analysis/data_merged/quenched/continuum_limit/logs/chi_dof.txt")

n = chi_dofs.shape[0]
av_chi_dof=0
min_chi_dof=99999999
max_chi_dof=0
for i in range(0,n):
    av_chi_dof += chi_dofs[i,1]
    if chi_dofs[i,1] < min_chi_dof:
        min_chi_dof = chi_dofs[i,1]
    if chi_dofs[i,1] > max_chi_dof:
        max_chi_dof = chi_dofs[i,1]
av_chi_dof /= 1.0*n
print(chi_dofs)
print("min     chi_dof: ", min_chi_dof) 
print("average chi_dof: ", av_chi_dof)
print("max     chi_dof: ", max_chi_dof) 
