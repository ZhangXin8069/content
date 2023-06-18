# get_ipython().run_line_magic('load_ext', 'blackcellmagic')

from iog_reader import iog_read
from os import listdir
from re import match
# import jax.numpy as np
import numpy as np
import pandas as pd
import gvar as gv
import lsqfit
# import matplotlib.pyplot as plt

pd.set_option("display.max_rows", None)

# @ti.kernel


# In[28]:


def get_c2pt(
    iog_path="Data",
    intrptr=["cnfg", "hdrn", "t", "snksmr", "lnk", "tsrc", "mntm"],
    N=24,
):
    iog_files =[
        iog_path + "/" + file
        for file in listdir(iog_path)
        if match(r".*dat.iog", file) != None
    ]
    df_list = []
    for iog_file in iog_files:
        df0 = iog_read(iog_file, intrptr)
        # print(
        #     df0[
        #         (df0["tsrc"] == "505046")
        #         & (df0["lnk"] == "9200505")
        #         & (df0["snksmr"] == "1")
        #     ][["Re", "cnfg", "hdrn", "t", "snksmr", "lnk", "tsrc", "mntm"]]
        # )
        df1 = df0[
            (df0["tsrc"] == "505046")
            & (df0["lnk"] == "9200505")
            & (df0["snksmr"] == "1")
        ]["Re"]
        df1 = df1.reset_index(drop=True)
        df2 = np.array([df1[i::N] for i in range(N)])
        df3 = np.mean(df2, axis=1)
        df_list.append(df3)
    print(np.array(df_list).shape)
    return np.array(df_list)


def jcknf_c2pt(c2pt=None):
    Ncnfg = c2pt.shape[0]
    c2pt_jcknf = (np.sum(c2pt, axis=0) - c2pt) / (Ncnfg - 1)
    return c2pt_jcknf


def make_mdls(t_dctnry, p):
    mdls = {}
    ts = t_dctnry["c2pt"]
    mdls["c2pt"] = p["n0"] * np.exp(p["n1"] + p["n2"] * ts)
    return mdls


def fit_c2tp(c2pt_jcknf, make_mdls):
    Ncnfg = c2pt_jcknf.shape[0]
    T = c2pt_jcknf.shape[1]
    T_start = 1
    T_end = T
    t_ary = np.array(range(T_start, T_end))
    c2pt_jcknf_avg = c2pt_jcknf[:, T_start:T_end]
    c2pt_avg_cntrl = np.mean(c2pt_jcknf_avg, axis=0)
    c2pt_avg_cov = (Ncnfg - 1) * np.cov(np.transpose(c2pt_jcknf_avg, axes=(1, 0)))
    tsep_dctnry = {"c2pt": t_ary}
    c2pt_dctnry = {"c2pt": gv.gvar(c2pt_avg_cntrl, c2pt_avg_cov)}

    p0 = {
        "n0": 1,
        "n1": -25,
        "n2": -3.2,
    }
    fit = lsqfit.nonlinear_fit(data=(tsep_dctnry, c2pt_dctnry), fcn=make_mdls, p0=p0)
    print(fit.format(True))


# In[6]:


#     t_ary = fit.data[0]["c2pt"]
#     c2pt_avg_dat_gvar = fit.data[1]["c2pt"]
#     c2pt_avg_dat_cntrl = np.array([c2.mean for c2 in c2pt_avg_dat_gvar])
#     c2pt_avg_dat_err = np.array([c2.sdev for c2 in c2pt_avg_dat_gvar])
#     t_lst = np.linspace(10, T - 10, 50)
#     c2pt_fit_fcn_gvar = fit.fcn({"c2pt": t_lst}, fit.p)["c2pt"]
#     c2pt_fit_fcn_cntrl = np.array([c2.mean for c2 in c2pt_fit_fcn_gvar])
#     c2pt_fit_fcn_err = np.array([c2.sdev for c2 in c2pt_fit_fcn_gvar])
#     c2pt_cntrl = np.mean(c2pt_jcknf_avg, axis=0)
#     c2pt_err = np.sqrt(Ncnfg - 1) * np.std(c2pt_jcknf_avg, axis=0)


#     plt.figure(dpi=600)
#     plt.errorbar(
#         t_ary,
#         c2pt_cntrl,
#         yerr=c2pt_err,
#         fmt="bo",
#         label="$C_2$",
#     )
#     plt.errorbar(
#         t_ary,
#         c2pt_avg_dat_cntrl,
#         yerr=c2pt_avg_dat_err,
#         fmt="go",
#         label="frwrd/bckwrd avg. $C_2$",
#     )
#     plt.plot(t_lst, c2pt_fit_fcn_cntrl, color="b", label="best fit")
#     plt.fill_between(
#         t_lst,
#         c2pt_fit_fcn_cntrl - c2pt_fit_fcn_err,
#         c2pt_fit_fcn_cntrl + c2pt_fit_fcn_err,
#     )
#     plt.xlabel("t/a")
#     plt.ylabel("$C_2$")
#     plt.legend(
#         loc="upper center",
#         frameon=True,
#         fancybox=True,
#         markerfirst=True,
#     )
#     plt.savefig("c2pt_fit.png")
#     plt.show()


# In[29]:


def main():
    c2pt=get_c2pt(iog_path='/public/home/sunpeng/chen_c/22_04_13_run_chroma/L24x72_m0235_3pt_mom4_st_20_sm_link/data')
    #c2pt = get_c2pt(iog_path="Data")
    c2pt_jcknf = jcknf_c2pt(c2pt=c2pt)
    fit_c2tp(c2pt_jcknf=c2pt_jcknf, make_mdls=make_mdls)


# get_ipython().run_line_magic('time', 'throws =main()')
if __name__ == "__main__":
    main()


# In[ ]:




from iog_reader import iog_read
