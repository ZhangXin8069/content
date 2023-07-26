#!/usr/bin/env python3
# encoding: utf-8
from iog_reader import iog_read

iog_file = "./Data/2pt_5350.dat.iog"
intrptr = ["cnfg","hdrn","t"]

iog = iog_read(iog_file,intrptr)
print(iog,'\n--------\n')

iog_pi=iog.loc[(iog['hdrn']=='0')]
print(iog_pi['Re'].to_numpy(),'\n--------\n')
print(iog_pi['Re'].to_list(),'\n--------\n')
