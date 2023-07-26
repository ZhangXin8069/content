import numpy as np
import math
import random as rd




def bootstrap(a,nconf): 
    conf,Tmax,parts = a.shape
    list = [rd.randint(0,conf-1) for i in range(nconf * conf)] 
    long = a[list,...]
    boot = np.reshape(long,(nconf,conf,Tmax,parts))
    return np.mean(boot,1) 

def jackknife_D3(a): 
    con,t,parts = a.shape
    ave = np.broadcast_to(np.mean(a,0),(con,t,parts)) 
    jac = (ave * con - a) / (con - 1) 
    return jac 

def jackknife_D4(a): 
    con,c,t,parts = a.shape
    ave = np.broadcast_to(np.mean(a,0),(con,c,t,parts)) 
    jac = (ave * con - a) / (con - 1) 
    return jac 
