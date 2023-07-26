#!/usr/bin/env python3
# encoding: utf-8


import numpy as np
from os import listdir
from iog_reader import iog_read
from re import match
import cProfile
from functools import reduce


iog_path = "./Data/"
T = 72

files = listdir(iog_path)
iogs = [file for file in files if match(r".*2pt.*.dat.iog",file)!=None]
Ncnfg = len(iogs)

c2pt = np.zeros((Ncnfg,T))
for indx in range(0,Ncnfg):
    dat=iog_read(iog_path+iogs[indx],["cnfg","hdrn","t"])
    c2pt[indx] = dat.loc[(dat["hdrn"]=="0")]["Re"].to_numpy()


c2pt_jcknf_1 = np.zeros((Ncnfg,T))
c2pt_jcknf_2 = np.zeros((Ncnfg,T))
c2pt_jcknf_3 = np.zeros((Ncnfg,T))

def c2pt_jcknf_loop_1():
    for i in range(0,Ncnfg):
        for j in range(0,T):
            c2pt_jcknf_1[i,j] = np.mean(np.delete(c2pt,[i],axis=0)[:,j])

def c2pt_jcknf_loop_2():
    for i in range(0,Ncnfg):
        c2pt_jcknf_2[i] = np.mean(np.delete(c2pt,[i],axis=0),axis=0)
    
def c2pt_jcknf():
    c2pt_sum = np.sum(c2pt,axis=0)
    global c2pt_jcknf_3 
    c2pt_jcknf_3 = ((c2pt_sum) - c2pt)/(Ncnfg-1)


cProfile.run('c2pt_jcknf_loop_1()')
print(" ------ \n")
cProfile.run('c2pt_jcknf_loop_2()')
print(" ------ \n")
cProfile.run('c2pt_jcknf()')
print(" ------ \n")

print(reduce(lambda x,y: x&y,np.fabs((c2pt_jcknf_1-c2pt_jcknf_2)/c2pt_jcknf_2)<1e-8))
print(reduce(lambda x,y: x&y,np.fabs((c2pt_jcknf_2-c2pt_jcknf_3)/c2pt_jcknf_2)<1e-8))
