#!/usr/bin/env python3
# encoding: utf-8
from ctypes import *
import os
import numpy as np
import pandas as pd

pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows', None)

lib = cdll.LoadLibrary(os.getcwd()+ '/iog.so')

nt = 72

def main():
    ## save 1file to txt
    ##filename = "../Data/2pt_4005.dat.iog"
    ##dataFrm = iog_read(filename)
    ##dataFrm.to_csv('data.txt', sep='\t')

    readin_data(tseq='3', save_file='yes')
    readin_data(tseq='4', save_file='yes')
    readin_data(tseq='5', save_file='yes')
    readin_data(tseq='6', save_file='yes')
    readin_data(tseq='7', save_file='yes')
    #eff_mass()
    

def iog_read(iog_fl):
    iog_file = c_char_p((iog_fl).encode('utf-8'))
    size = lib.getsize(iog_file)

    class iog_type(Structure):
        _fields_ = [('info_shape', POINTER(c_int)), ('info', POINTER(c_int)), ('Re', POINTER(c_double)),  ('Im', POINTER(c_double))]


    lib.getdat.restype = POINTER(iog_type)
    iog = lib.getdat(iog_file)

    info_shape = np.ctypeslib.as_array(iog.contents.info_shape, shape=[2])
    info = np.ctypeslib.as_array(iog.contents.info, shape=[np.prod(info_shape)])
    re = np.ctypeslib.as_array(iog.contents.Re, shape=[size])
    im = np.ctypeslib.as_array(iog.contents.Im, shape=[size])
    lables = info.reshape(info_shape)

    iog_dat = np.concatenate((re.reshape(size,1), im.reshape(size,1)),axis=1)
    lables_str = np.vectorize(str)(lables)
    intrprtrs = ['cfg', 'hadrons', 't']

    iog_dict = {}
    for indx in range(0, len(intrprtrs)):
        iog_dict.update({intrprtrs[indx]:lables_str[:,indx]})

    iog_dict.update({'Re':re})
    iog_dict.update({'Im':im})
    #print(iog_dict)

    dataFrm = pd.DataFrame(iog_dict)
    return dataFrm

def readin_data(tseq, save_file='no'):
    cfg_list = []
    with open('../list_all', 'r') as f:
        for line in f.readlines():
            cfg_list.append(line.strip().split())
    cfg_list = np.squeeze(np.array(cfg_list, int)) # shape=(n_cfg, )

    data_collect = []
    for cfg in cfg_list:
        filename = "../Data/3pt_"+str(cfg)+"_tseq"+tseq+".dat.iog"
        dataFrm = iog_read(filename)
        data = np.array(dataFrm.values.tolist(), float)
        
        col_no = data.shape[-1] # number of columns
        data_reshape = data.reshape(-1, nt, col_no)
        data_collect.append(data_reshape)
    data_collect = np.array(data_collect) # shape=(n_cfg, n_hadron, nt, col_no)

    if(save_file=='yes'):
        np.save('data_3pt_tseq'+tseq+'.npy', data_collect)

    return data_collect
        
def eff_mass():
    data = np.load('data_2pt.npy')
    print(data.shape)

if __name__ == "__main__":
    main()


