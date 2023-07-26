#!/usr/bin/env python3
# encoding: utf-8

import numpy as np
import gvar as gv
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

alttc = 0.108
fm2GeV = 0.197

c2pt = np.load("./c2pt.npy")

Ncnfg = c2pt.shape[0]
T = c2pt.shape[1]
T_hlf = T//2


c2pt_sum = np.sum(c2pt,axis=0)
c2pt_jcknf = ((c2pt_sum) - c2pt)/(Ncnfg-1) # jack-knife resample
c2pt_cntrl = np.mean(c2pt_jcknf,axis=0) # jack-knife mean
c2pt_cov = (Ncnfg-1)*np.cov(np.transpose(c2pt_jcknf,axes=(1,0))) # jack-knife covariance
c2pt_err = np.sqrt(Ncnfg-1)*np.std(c2pt_jcknf,axis=0)



fig, ax = plt.subplots(1,1, figsize=(10, 6.18))
ax.errorbar(np.array(range(0,T)), c2pt_cntrl, yerr=c2pt_err, fmt = 'bo')
ax.set_yscale('log')
fig.savefig("./c2pt_dat.png")

c2pt_fwrd_jcknf = c2pt_jcknf[:,0:T_hlf+1]
c2pt_bwrd_jcknf = c2pt_jcknf[:,T_hlf:][:,::-1]
c2pt_avg_jcknf = (c2pt_fwrd_jcknf[:,1:]+c2pt_bwrd_jcknf)/2


m_eff_log_fwrd_jcknf = np.log(c2pt_fwrd_jcknf[:,:-1]/c2pt_fwrd_jcknf[:,1:])*fm2GeV/alttc
m_eff_log_bwrd_jcknf = np.log(c2pt_bwrd_jcknf[:,:-1]/c2pt_bwrd_jcknf[:,1:])*fm2GeV/alttc


m_eff_log_fwrd_cntrl = np.mean(m_eff_log_fwrd_jcknf,axis=0)
m_eff_log_fwrd_err = np.sqrt(Ncnfg-1)*np.std(m_eff_log_fwrd_jcknf,axis=0)
m_eff_log_bwrd_cntrl = np.mean(m_eff_log_bwrd_jcknf,axis=0)
m_eff_log_bwrd_err = np.sqrt(Ncnfg-1)*np.std(m_eff_log_bwrd_jcknf,axis=0)

m_eff_log_cntrl= np.concatenate((m_eff_log_fwrd_cntrl, m_eff_log_bwrd_cntrl[::-1]), axis=0)
m_eff_log_err= np.concatenate((m_eff_log_fwrd_err, m_eff_log_bwrd_err[::-1]), axis=0)

t_ary = np.array(range(0,T-1))
ini = 0.2*np.ones_like(t_ary)

def eff_mass_eqn(c2pt):
    return lambda E0: (c2pt[:-1]/c2pt[1:]) - np.cosh(E0*(T/2-t_ary)) / np.cosh(E0*(T/2-(t_ary+1)))
def fndroot(eqnf,ini):
    sol = fsolve(eqnf,ini, xtol=1e-5)
    return sol

m_eff_cosh_jcknf = np.array([fndroot(eff_mass_eqn(c2pt),ini) for c2pt in c2pt_jcknf])*fm2GeV/alttc
m_eff_cosh_cntrl = np.mean(m_eff_cosh_jcknf,axis=0)
m_eff_cosh_err = np.sqrt(Ncnfg-1)*np.std(m_eff_cosh_jcknf,axis=0)

xshft = 0.3
fig, ax = plt.subplots(1,1, figsize=(10, 7*0.5))
ax.set_ylim([-0.05, 0.4])
fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
ax.errorbar(t_ary+xshft, m_eff_log_cntrl, yerr=m_eff_log_err, fmt = 'go',alpha=0.5)
ax.errorbar(t_ary, m_eff_cosh_cntrl, yerr=m_eff_cosh_err, fmt = 'bo',alpha=0.5)
ax.set_xlabel('t/a')
ax.set_ylabel('$m_{\mathrm{eff}}$ [GeV]')
fig.savefig("./m_eff.png")
