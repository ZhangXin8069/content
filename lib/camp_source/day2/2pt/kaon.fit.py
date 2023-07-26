import numpy as np
import gvar as gv
from lsqfit import nonlinear_fit as nlinefit
import scipy.optimize as opt
import matplotlib.pyplot as plt
import seaborn as sns

sns.set()  # Set style for matplotlib

# Load data from .npy file
c2pt_list = np.load('./c2pt.npy')
c2pt_list = np.abs(c2pt_list)
# Get number of c2pt and length of time
NUMBER, TLENTH = c2pt_list.shape


def jackKnife(data: np.ndarray, dropfirst: bool, average: bool, normalize: bool) -> np.ndarray:
    """
    Perform jackknife resampling on input data.

    Args:
    - data: Input data of shape (n, m), where n is the number of samples and m is the length of time.
    - dropfirst: Whether to drop the first time slice of the data.
    - average: Whether to average over the time slices of the data.
    - normalize: Whether to normalize the data by dividing by the mean of the first time slice.

    Returns:
    - jackkniff_data: Jackknifed data of shape (n, m).
    - mean_jackkniff_data: Mean of jackknifed data over samples.
    - jackkniff_data_error: Error of jackknifed data over samples.
    - jackkniff_data_covar: Covariance matrix of jackknifed data over samples.
    """
    number = data.shape[0]
    if dropfirst:
        time_array = [i for i in range(1, data.shape[1])]
    else:
        time_array = [i for i in range(0, data.shape[1])]
    data = data[:, time_array]
    if average:
        time_half = int((data.shape[1] + 1) // 2)
        data = (data[:, :time_half] +
                data[:, data.shape[1] - time_half:][:, ::-1]) / 2.0
    if normalize:
        data /= np.mean(data[:, 0], axis=0)
    jackkniff_data = (np.sum(data, axis=0) - data) / (number - 1)
    mean_jackkniff_data = np.mean(jackkniff_data, axis=0)
    jackkniff_data_error = np.sqrt(number - 1) * np.std(jackkniff_data, axis=0)
    jackkniff_data_covar = (number - 1) * \
        np.cov(np.transpose(jackkniff_data, axes=(1, 0)))
    return jackkniff_data, mean_jackkniff_data, jackkniff_data_error, jackkniff_data_covar


# Preprocess data using jackknife resampling
jackkniff_c2pt_list, mean_jackkniff_c2pt_list, jackkniff_c2pt_list_error, jackkniff_c2pt_list_covar = jackKnife(
    data=c2pt_list, dropfirst=True, average=True, normalize=True)


def modelFunc(t: np.ndarray, p: dict) -> np.ndarray:
    """
    Define the model function for fitting.

    Args:
    - t: Time array of shape (m,) for fitting.
    - p: Parameters dictionary with keys 'C0' and 'E0'.

    Returns:
    - Model function evaluated at t with given parameters.
    """
    return p['C0'] * np.cosh(p['E0'] * (t - TLENTH / 2.0))


def nonlinearFit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray, modelfunc: callable) -> list:
    """
    Perform nonlinear fitting on input data.

    Args:
    - x: Independent variable array of shape (m,) for fitting.
    - y: Dependent variable array of shape (m,) for fitting.
    - yerr: Error array of shape (m, m) for fitting.
    - modelfunc: Model function to fit the data.

    Returns:
    - params: Fitted parameter values.
    - errors: Fitted parameter errors.
    - chi2: Reduced chi-squared value.
    - logGBF: Logarithm of the generalized Bayesian factor.
    """
    fit = nlinefit(data=(x, gv.gvar(y, yerr)), prior=gv.gvar(
        dict(C0='1(100)', E0='1(100)')), fcn=modelfunc, debug=True)
    print(fit.format(True))
    params = [i.val for i in fit.p.values()]
    errors = [i.sdev for i in fit.p.values()]
    chi2 = fit.chi2 / fit.dof
    logGBF = fit.logGBF
    return params, errors, chi2, logGBF

print("*********KAON*********")
start = 9
end = 24

params, errors, chi2, logGBF = nonlinearFit(x=np.array(range(start, end)),
                                            y=mean_jackkniff_c2pt_list[start:end],
                                            yerr=jackkniff_c2pt_list_covar[start:end,
                                                                           start:end],
                                            modelfunc=modelFunc)
print(params,chi2)
x = np.array(range(1, TLENTH))
_, y, yerr, _ = jackKnife(
    data=c2pt_list, dropfirst=False, average=False, normalize=True)
plt.errorbar(x=x, y=y[x], yerr=yerr[x], label='$C_2$')
plt.fill_between(x=x, y1=y[x]-yerr[x], y2=y[x]+yerr[x], alpha=0.3)
plt.plot(x, modelFunc(x, dict(C0=params[0], E0=params[1])), label='$best fit$')
plt.xlabel('t/a')
plt.ylabel('$C_2$')
plt.yscale("log")
plt.legend()
plt.title("KAON FIT")
plt.savefig("./kaon.fit.png")
