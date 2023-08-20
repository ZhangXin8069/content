import numpy as np
import gvar as gv
from lsqfit import nonlinear_fit as nlinefit
import matplotlib.pyplot as plt
import seaborn as sns
sns.set()

# init parameters
np.random.seed(101)
NUMBER = 10000
TLENTH = 100
N_DUMP = 1000
N_SKIP = 20
DELTA = 0.05
A = 0.1
data_list = np.zeros((NUMBER, TLENTH))  # cold start

# define some functions


def S(data: np.ndarray, t: int):
    b_t = (t-1+TLENTH) % TLENTH
    f_t = (t+1) % TLENTH
    return (data[t]-data[b_t])**2.0/(2.0*A)+A*(data[t]**2+data[b_t]**2)/4.0+(data[f_t]-data[t])**2.0/(2.0*A)+A*(data[f_t]**2+data[t]**2)/4.0


def jackKnife(data_list: np.ndarray) -> np.ndarray:
    number = data_list.shape[0]
    jackkniff_data_list = (np.sum(data_list, axis=0) -
                           data_list) / (number - 1)
    mean_jackkniff_data_list = np.mean(jackkniff_data_list, axis=0)
    jackkniff_data_list_error = np.sqrt(
        number - 1) * np.std(jackkniff_data_list, axis=0)
    jackkniff_data_list_covar = (number - 1) * \
        np.cov(np.transpose(jackkniff_data_list, axes=(1, 0)))
    return jackkniff_data_list, mean_jackkniff_data_list, jackkniff_data_list_error, jackkniff_data_list_covar


def modelFunc(t: np.ndarray, p: dict) -> np.ndarray:
    return p['C0'] * np.exp(p['E0'] * t)


def nonlinearFit(x: np.ndarray, y: np.ndarray, yerr: np.ndarray, modelfunc: callable) -> list:
    fit = nlinefit(data=(x, gv.gvar(y, yerr)), prior=gv.gvar(
        dict(C0='1(100)', E0='1(100)')), fcn=modelfunc, debug=True)
    print(fit.format(True))
    params = [i.val for i in fit.p.values()]
    errors = [i.sdev for i in fit.p.values()]
    chi2 = fit.chi2 / fit.dof
    logGBF = fit.logGBF
    return params, errors, chi2, logGBF


# generating configurations
data = np.zeros(TLENTH)
for i in range(NUMBER):
    t_list = np.arange(TLENTH)
    # np.random.shuffle(t_list)
    for t in t_list:
        data_new = data+np.random.uniform(-1, 1)*DELTA
        if (np.random.uniform(0, 1) < np.exp(S(data=data, t=t)-S(data=data_new, t=t))):
            data[t] = data_new[t]
        else:
            data_new[t] = data[t]
    data_list[i, :] = data_new
data_list = data_list[N_DUMP:NUMBER:N_SKIP, :]

# measuring the ground state energy
data2_list = data_list**2
jackkniff_data2_list, mean_jackkniff_data2_list, jackkniff_data2_list_error, jackkniff_data2_list_covar = jackKnife(
    data_list=data2_list)
ground_energy = np.mean(jackkniff_data2_list, axis=1)
ground_energy_mean = np.mean(ground_energy)
ground_energy_error = np.std(ground_energy)*(ground_energy.shape[0]-1)**0.5

# give Theoretical value
theo_data = np.zeros(TLENTH)
R = 1+A**2/4.0-A*(1+A**2/4.0)**0.5
C = R*0.5*(1+A**2/4.0)**0.5
for t in range(TLENTH):
    try:
        theo_data[t] = C*(1+R**t)/(1-R**t)
    except ZeroDivisionError:
        theo_data[t] = C*0.5

print("*********HARM_OSCIL_1D_GROUND_ENERGY*********")
x = np.array(range(TLENTH))
plt.plot(x, mean_jackkniff_data2_list,
         label="mean_jackkniff_data2_list", linewidth=1, color='blue')
plt.plot(x[3:], theo_data[x[3:]], label="theo_data", linewidth=4, color='red')
plt.errorbar(
    x=x, y=mean_jackkniff_data2_list, yerr=jackkniff_data2_list_error)
plt.fill_between(x=x, y1=mean_jackkniff_data2_list-jackkniff_data2_list_error,
                 y2=mean_jackkniff_data2_list+jackkniff_data2_list_error, alpha=0.3)
plt.xlabel("n")
plt.ylabel("<$u^2$>")
plt.legend()
plt.title("HARM_OSCIL_1D_GROUND_ENERGY\n"+"ground_energy_mean:" +
          str(ground_energy_mean)+"\nground_energy_error:"+str(ground_energy_error))
plt.savefig("HARM_OSCIL_1D_GROUND_ENERGY.png")
plt.show()

# measuring the excited state energy

c_data_list = []
data_list_top = data_list[:, :13]
for n in range(data_list_top.shape[0]):
    data_tmp = []
    for t in range(data_list_top.shape[1]):
        tmp = 0
        for t0 in range(data_list_top.shape[1]):
            tmp += data_list_top[n, (t+t0) %
                                 data_list_top.shape[1]]*data_list_top[n, t0]
            print(tmp)
        data_tmp.append(tmp/data_list_top.shape[1])
    c_data_list.append(data_tmp)
c_data_list = np.array(c_data_list)
jackkniff_c_data_list, mean_jackkniff_c_data_list, jackkniff_c_data_list_error, jackkniff_c_data_list_covar = jackKnife(
    data_list=c_data_list)
x = np.array(range(jackkniff_c_data_list.shape[1]))
params, errors, chi2, logGBF = nonlinearFit(
    x=x, y=mean_jackkniff_c_data_list, yerr=jackkniff_c_data_list_error, modelfunc=modelFunc)
print(mean_jackkniff_c_data_list)
print("*********HARM_OSCIL_1D_EXCITED_STATE_ENERGY*********")
plt.plot(x, mean_jackkniff_c_data_list,
         label="mean_jackkniff_c_data_list", linewidth=1, color='blue')
plt.errorbar(
    x=x, y=mean_jackkniff_c_data_list, yerr=jackkniff_c_data_list_error)
plt.fill_between(x=x, y1=mean_jackkniff_c_data_list-jackkniff_c_data_list_error,
                 y2=mean_jackkniff_c_data_list+jackkniff_c_data_list_error, alpha=0.3)
plt.xlabel("t")
plt.ylabel("<$C(t)$>")
plt.yscale("log")
plt.legend()
# plt.title("HARM_OSCIL_1D_GROUND_ENERGY\n"+"ground_energy_mean:" +
#           str(ground_energy_mean)+"\nground_energy_error:"+str(ground_energy_error))
plt.savefig("HARM_OSCIL_1D_EXCITED_STATE_ENERGY.png")
plt.show()