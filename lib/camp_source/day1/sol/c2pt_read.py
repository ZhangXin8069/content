#!/usr/bin/env python3
# encoding: utf-8

import numpy as np
from os import listdir
from iog_reader import iog_read
from re import match
import cProfile
from functools import reduce


iog_path = "./Data/"

files = listdir(iog_path)
iogs = [file for file in files if match(r".*2pt.*.dat.iog",file)!=None]
Ncnfg = len(iogs)
T = 72

c2pt = np.zeros((Ncnfg,T))
for indx in range(0,Ncnfg):
    dat=iog_read(iog_path+iogs[indx],["cnfg","hdrn","t"])
    c2pt[indx] = dat.loc[(dat["hdrn"]=="0")]["Re"].to_numpy()


np.save("c2pt.npy",c2pt)