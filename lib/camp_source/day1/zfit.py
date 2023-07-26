import numpy as np
import gvar as gv
from lsqfit import nonlinear_fit as nlinefit
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()  # beautiful

c2pt_list = np.load("./c2pt.npy")  # load from .npy
NUMBER, TLENTH = c2pt_list.shape  # number of c2pt,lenth of time


def jackKnife(data: np.ndarray, dropfirst: bool, average: bool, normalize: bool) -> np.ndarray:
    number = data.shape[0]
    if(dropfirst):
        time_array = [i for i in range(1, data.shape[1]+1)]
    else:
        time_array = [i for i in range(0, data.shape[1]+1)]
    data = data[:, time_array]
    if(average):
        time_half = int((data.shape[1]+1)//2)
        data = (data[:, :time_half] +
                data[:, data.shape[1] - time_half:][:, ::-1])/2.0
    if(normalize):
        data /= np.mean(data[:, 0], axis=0)
    jackkniff_data = (np.sum(data, axis=0) - data) / (number - 1)
    mean_jackkniff_data = np.mean(jackkniff_data, axis=0)
    jackkniff_data_error = np.sqrt(number-1)*np.std(jackkniff_data, axis=0)
    jackkniff_data_covar = (number-1) * \
        np.cov(np.transpose(jackkniff_data, axes=(1, 0)))
    return jackkniff_data, mean_jackkniff_data, jackkniff_data_error, jackkniff_data_covar


# data preprocessing
jackkniff_c2pt_list, mean_jackkniff_c2pt_list, jackkniff_c2pt_list_error, jackkniff_c2pt_list_covar = jackKnife(
    data=c2pt_list, dropfirst=True, average=True, normalize=True)


def modelFunc(t, p):
    return p["C0"] * np.cosh(p["E0"] * (t-TLENTH/2.0))

def nonlinearFit(x, y, yerr, modelfunc):
    prior = gv.gvar(dict(C0='0(5)', E0='0.2(0.5)'))
    fit = nlinefit(
        data=(x, gv.gvar(y, yerr)), prior=prior, fcn=modelfunc, debug=True)
    print(fit.format(True))
    return [i.val for i in fit.p.values()], [i.sdev for i in fit.p.values()], fit.chi2/fit.dof, fit.logGBF

T_strt = 10 #? starting point of the fit
T_end =  30 #? ending point of the fit
x=[i for i in range(T_strt, T_end+1)]
tsep_dctnry = {'c2pt': t_ary[T_strt:T_end]}
c2pt_dctnry = {'c2pt': gv.gvar(c2pt_avg_cntrl[T_strt:T_end], c2pt_avg_cov[T_strt:T_end, T_strt:T_end])