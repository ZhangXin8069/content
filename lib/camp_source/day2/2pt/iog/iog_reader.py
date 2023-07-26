#!/usr/bin/env python3
# enCoding: utf-8

#?read iog file into pandas dataframe

from ctypes import *
import os
import numpy as np
import pandas as pd


lib = cdll.LoadLibrary(os.getcwd()+ '/iog/iog.so')


def iog_read(iog_fl, intrprtrs):
    iog_file = c_char_p(iog_fl.encode('utf-8'))
    
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

    lables_str = np.vectorize(str)(lables)

    iog_dict = {}
    for indx in range(0, len(intrprtrs)):
        iog_dict.update({intrprtrs[indx]:lables_str[:,indx]})

    iog_dict.update({'Re':re})
    iog_dict.update({'Im':im})

    dataFrm = pd.DataFrame(iog_dict)
    return dataFrm
