import numpy as np
import matplotlib.pyplot as plt
import gvar as gv
from lsqfit import nonlinear_fit as nlinefit

# Init parameters
N_conf = 100000
N_dump = 50000
N_skip = 100
M = int((N_conf-N_dump)/N_skip)
T = 50
N_t = 500
a = T/N_t
delta = 0.05

# Define some functions


def jaknf(u):
    return (np.sum(u, axis=0)-u)/(u.shape[0]-1)


def modelFunc(t, p: dict):
    return p['C0'] * np.exp(-p['E0'] * t)


def nonlinearFit(x, y, yerr, modelfunc: callable):
    fit = nlinefit(data=(x, gv.gvar(y, yerr)), prior=gv.gvar(
        dict(C0='1(100)', E0='1(100)')), fcn=modelfunc, debug=True)
    print(fit.format(True))
    params = [i.val for i in fit.p.values()]
    errors = [i.sdev for i in fit.p.values()]
    chi2 = fit.chi2 / fit.dof
    logGBF = fit.logGBF
    return params, errors, chi2, logGBF


def weight(u, t):
    u_b_t = u[(t-1+N_t) % N_t]
    u_f_t = u[(t+1) % N_t]
    u_t = u[t]
    S = (u_f_t-u_t)**2/(2.0*a)+a*(u_f_t**2+u_t**2)/4.0 + \
        (u_t-u_b_t)**2/(2.0*a)+a*(u_t**2+u_b_t**2)/4.0
    return np.exp(-1.0*S)


# Generating configurations
Un = np.zeros((N_conf, N_t))
old_u = np.zeros(N_t)
new_u = np.zeros(N_t)


for n in range(N_conf):
    t_iter = np.array(range(N_t))
    np.random.shuffle(t_iter)
    for t in t_iter:
        new_u = old_u+delta*np.random.uniform(-1, 1)
        if np.random.uniform(0, 1) < weight(new_u, t)/weight(old_u, t):
            old_u[t] = new_u[t]
        else:
            new_u[t] = old_u[t]
    Un[n, :] = new_u
Um = Un[N_dump:N_conf:N_skip, :]
# Um_2 = Um**2
# jk_Um_2 = jaknf(Um_2)

# u_2 = np.mean(np.mean(jk_Um_2, axis=1), axis=0)
# u_2_err = np.sum((np.mean(jk_Um_2, axis=1)-u_2)**2)**0.5
# print("theory u_2:"+str(0.5*(1+a**2/4)**0.5) +
#       "\nu_2:"+str(u_2)+"\nu_2_err:" + str(u_2_err))
# theory u_2:0.5006246098625197
# u_2:0.5013973234527527
# u_2_err:0.01606474433280273
np.save("Um.npy", Um)
# Um = np.load("Um.npy")

# Give Theoretical value
theo_u_2 = np.zeros(N_t)
R = 1+a**2/4.0-a*(1+a**2/4.0)**0.5
C = 0.5*(1+a**2/4.0)**0.5
for t in range(N_t):
    try:
        theo_u_2[t] = C*(1+R**t)/(1-R**t)
    except ZeroDivisionError:
        theo_u_2[t] = C
# Measuring the ground energy
print("*********HARM_OSCIL_1D_GROUND_ENERGY*********")
t_iter = np.array(range(N_t))
Um_2 = Um**2
jk_Um_2 = jaknf(Um_2)
m_iter = np.array(range(N_conf))[N_dump:N_conf:N_skip]
u_2_m = np.mean(jk_Um_2, axis=1)
plt.scatter(m_iter, u_2_m)
# jk_Um_2_mean = np.mean(jk_Um_2, axis=0)
# jk_Um_2_err = (M-1)**0.5*np.std(jk_Um_2, axis=0)
# u_2 = np.mean(np.mean(jk_Um_2, axis=1), axis=0)
# u_2_err = np.sum((np.mean(jk_Um_2, axis=1)-u_2)**2)**0.5
# plt.scatter(t_iter, jk_Um_2_mean,
#             label="jk_Um_2_mean", color='blue')
# plt.scatter(t_iter[3:], theo_u_2[t_iter[3:]],
#             label="theo_u_2", color='red')
# plt.errorbar(
#     x=t_iter, y=jk_Um_2_mean, yerr=jk_Um_2_err)
# plt.fill_between(x=t_iter, y1=jk_Um_2_mean-jk_Um_2_err,
#                  y2=jk_Um_2_mean+jk_Um_2_err, alpha=0.3)
plt.xlabel("m")
plt.ylabel("<$u^2$>")
plt.legend()
# plt.title("HARM_OSCIL_1D_GROUND_ENERGY\n"+"ground_energy_mean:" +
#           str(u_2)+"\nground_energy_error:"+str(u_2_err))
plt.savefig("HARM_OSCIL_1D_GROUND_ENERGY.png")
plt.show()

# Measuring the excited state energy

# Ct = []
# jk_Um = jaknf(Um)
# for t in range(N_t):
#     tmp_list = []
#     for m in range(M):
#         tmp_list.append(
#             np.mean([Um[m, (t+t0) % N_t]*Um[m, t0] for t0 in range(N_t)]))
#     Ct.append(np.mean(tmp_list))
# Ct = np.array(Ct)
# np.save("Ct.npy", Ct)
# print(Ct)
# t_iter = np.array(range(N_t))
# Ct = np.load("Ct.npy")
# jk_Ct = jaknf(Ct)
# jk_Ct_mean = np.mean(jk_Ct, axis=0)
# jk_Ct_err = (M-1)**0.5*np.std(jk_Ct, axis=0)
# # params, errors, chi2, logGBF = nonlinearFit(
# #     x=t_iter, y=jk_Ct_mean, yerr=jk_Ct_err, modelfunc=modelFunc)
# # # print(jk_Ct_mean)
# print("*********HARM_OSCIL_1D_EXCITED_STATE_ENERGY*********")
# plt.plot(t_iter, jk_Ct_mean, label="jk_Ct_mean", linewidth=1, color='blue')
# plt.errorbar(
#     x=t_iter, y=jk_Ct_mean, yerr=jk_Ct_err)
# plt.fill_between(x=t_iter, y1=jk_Ct_mean-jk_Ct_err,
#                  y2=jk_Ct_mean+jk_Ct_err, alpha=0.3)
# plt.xlabel("t")
# plt.ylabel("<$C(t)$>")
# plt.yscale("log")
# plt.legend()
# plt.title("HARM_OSCIL_1D_EXCITED_STATE_ENERGY")
# plt.savefig("HARM_OSCIL_1D_EXCITED_STATE_ENERGY.png")
# plt.show()
