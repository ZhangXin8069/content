#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# %matplotlib inline
# %load_ext blackcellmagic
# %load_ext watermark
# %watermark
# %watermark -p openmc
from IPython.display import Image
import numpy as np
import openmc
from math import pi
import matplotlib.pyplot as plt
from matplotlib import pyplot
import openmc.lib
import openmc.deplete
eV=1.60217662e-19
path="/home/zhangxin/endfb71_hdf5/chain_simple.xml"


# In[ ]:


# x=0.198#keff=0.98*
# x=0.193#keff=0.96*
# x=0.192#keff=0.94*
import sys
x= float(sys.argv[1])


# In[ ]:


Tem0 = 20 + 273.15

nam_ele = ["Th_232", "U_233", "O", "Fe", "Cr", "Mn", "W", "Pb"]
num_ele = [3, 2, 3, 4, 4, 4, 4, 5]


def ele_nam(nam="O", num=4):
    return [nam + "_" + str(i) for i in range(1, num + 1)]


ele_nam_lis0 = [ele_nam(i, j) for i, j in zip(nam_ele, num_ele)]
ele_nam_lis = []
for j in ele_nam_lis0:
    for i in j:
        ele_nam_lis.append(i)
num_den = [
    6.35e-3 * (1 - x),
    7.45e-3 * (1 - x),
    7.45e-3,
    6.35e-3 * x,
    7.45e-3 * x,
    1.27e-2,
    1.49e-2,
    1.49e-2,
    8.10e-3,
    8.87e-3,
    8.87e-3,
    6.63e-3,
    1.12e-3,
    1.06e-3,
    1.06e-3,
    8.00e-4,
    4.60e-5,
    5.10e-5,
    5.10e-5,
    3.80e-5,
    4.60e-5,
    5.10e-5,
    5.10e-5,
    3.80e-5,
    1.77e-2,
    1.56e-2,
    1.56e-2,
    2.41e-2,
    3.05e-2,
]

ele_set = {i: j for i, j in zip(ele_nam_lis, num_den)}


# In[ ]:


mat_lis = [0, 1, 2, 3, 4]
for i in mat_lis:
    mat_lis[i] = openmc.Material(temperature=Tem0)
    lis = []
    for j in ele_nam_lis:
        if j.split(sep="_")[-1] == str(i + 1):
            lis.append(j)
    atom_sum = sum([ele_set[eva] for eva in lis])
    for k in lis:
        if len(k.split(sep="_")) == 2:
            mat_lis[i].add_element(
                element=str(k.split(sep="_")[0]),
                percent=ele_set[k] / atom_sum,
                percent_type="ao",
            )
        if len(k.split(sep="_")) == 3:
            mat_lis[i].add_nuclide(
                nuclide=str(k.split(sep="_")[0]) + str(k.split(sep="_")[1]),
                percent=ele_set[k] / atom_sum,
                percent_type="ao",
            )
    mat_lis[i].set_density(units="atom/cm3", density=atom_sum * 1e24)


# In[ ]:


mats_lis = [
    mat_lis[0],  # 1
    mat_lis[1],  # 2
    mat_lis[2],  # 3
    mat_lis[3],  # 4up
    mat_lis[3],  # 4do
    mat_lis[4],  # 5in
    mat_lis[4],  # 5mup
    mat_lis[4],  # 5mdo
    mat_lis[4],  # 5ou
]
mats = openmc.Materials(mat_lis)
rad = [0, 32.5, 82.5, 147.5, 162.5, 320]
hig = [-320, -165, -75, 75, 165, 320]
# ins out ups dow


# In[ ]:


ins = [
    openmc.ZCylinder(r=rad[1]),
    openmc.ZCylinder(r=rad[2]),
    openmc.ZCylinder(r=rad[3]),
    openmc.ZCylinder(r=rad[1]),
    openmc.ZCylinder(r=rad[1]),
    openmc.ZCylinder(r=rad[0]),
    openmc.ZCylinder(r=rad[1]),
    openmc.ZCylinder(r=rad[1]),
    openmc.ZCylinder(r=rad[4]),
]
out = [
    openmc.ZCylinder(r=rad[2]),
    openmc.ZCylinder(r=rad[3]),
    openmc.ZCylinder(r=rad[4]),
    openmc.ZCylinder(r=rad[4]),
    openmc.ZCylinder(r=rad[4]),
    openmc.ZCylinder(r=rad[1]),
    openmc.ZCylinder(r=rad[4]),
    openmc.ZCylinder(r=rad[4]),
    openmc.ZCylinder(r=rad[5]),
]
dow = [
    openmc.ZPlane(z0=hig[2]),
    openmc.ZPlane(z0=hig[2]),
    openmc.ZPlane(z0=hig[2]),
    openmc.ZPlane(z0=hig[3]),
    openmc.ZPlane(z0=hig[1]),
    openmc.ZPlane(z0=hig[0]),
    openmc.ZPlane(z0=hig[4]),
    openmc.ZPlane(z0=hig[0]),
    openmc.ZPlane(z0=hig[0]),
]
ups = [
    openmc.ZPlane(z0=hig[3]),
    openmc.ZPlane(z0=hig[3]),
    openmc.ZPlane(z0=hig[3]),
    openmc.ZPlane(z0=hig[4]),
    openmc.ZPlane(z0=hig[2]),
    openmc.ZPlane(z0=hig[5]),
    openmc.ZPlane(z0=hig[5]),
    openmc.ZPlane(z0=hig[1]),
    openmc.ZPlane(z0=hig[5]),
]


# In[ ]:


for i in ins + out + dow + ups:
    i.boundary_type = "transmission"
out[8].boundary_type = "vacuum"
dow[5].boundary_type = "vacuum"
dow[7].boundary_type = "vacuum"
dow[8].boundary_type = "vacuum"
ups[5].boundary_type = "vacuum"
ups[6].boundary_type = "vacuum"
ups[8].boundary_type = "vacuum"


# In[ ]:


ins_lis = [
    rad[1],
    rad[2],
    rad[3],
    rad[1],
    rad[1],
    rad[0],
    rad[1],
    rad[1],
    rad[4],
]
out_lis = [
    rad[2],
    rad[3],
    rad[4],
    rad[4],
    rad[4],
    rad[1],
    rad[4],
    rad[4],
    rad[5],
]
dow_lis = [
    hig[2],
    hig[2],
    hig[2],
    hig[3],
    hig[1],
    hig[0],
    hig[4],
    hig[0],
    hig[0],
]
ups_lis = [
    hig[3],
    hig[3],
    hig[3],
    hig[4],
    hig[2],
    hig[5],
    hig[5],
    hig[1],
    hig[5],
]


# In[ ]:


cel_lis = [0, 1, 2, 3, 4, 5, 6, 7, 8]
for i in cel_lis:
    mats_lis[i].volume = (
        (ups_lis[i] - dow_lis[i]) * (out_lis[i] ** 2 - ins_lis[i] ** 2) * pi
    )
    cel_lis[i] = openmc.Cell()
    cel_lis[i].region = +ins[i] & -out[i] & +dow[i] & -ups[i]
    cel_lis[i].fill = mats_lis[i]
    cel_lis[i].temperature = Tem0


# In[ ]:


imp = 500
cel_voi = openmc.Cell()
voi_spa = [
    +openmc.ZPlane(z0=hig[0] - imp),
    -openmc.ZPlane(z0=hig[5] + imp),
    -openmc.ZCylinder(r=rad[5] + imp),
]
cel_voi_space = voi_spa[0]
for i in voi_spa:
    i.boundary_type = "vacuum"
    cel_voi_space = cel_voi_space & i
cel_voi.region = cel_voi_space
cel_voi.temperature = Tem0
cel_voi.fill = mat_lis[4]
cel_lis.append(cel_voi)


# In[ ]:


uni = openmc.Universe(cells=cel_lis)
# root_cell = openmc.Cell(
#     fill=uni,
#     region=
#     +openmc.ZPlane(z0=hig[0])
#     & -openmc.ZPlane(z0=hig[5])
#     & -openmc.ZCylinder(r=rad[5]),
# )


# In[ ]:


# width=(640,640)
# pixels=(640,640)
# color_by='material'


# In[ ]:


# uni.plot(origin=(0  , 0, -1), width=width, basis="xy", pixels=pixels, color_by=color_by)


# In[ ]:


# uni.plot(origin=(0, 0, -76), width=width, basis="xy", pixels=pixels, color_by=color_by)


# In[ ]:


# uni.plot(origin=(0, 0,- 166), width=width, basis="xy", pixels=pixels, color_by=color_by)


# In[ ]:


# uni.plot(origin=(0, 1, 0), width=width, basis="xz", pixels=pixels, color_by=color_by)


# In[ ]:


geo = openmc.Geometry(uni)
# geo=openmc.Geometry([root_cell])


# In[ ]:


sou = openmc.Source()
sou.space = openmc.stats.Box(
    (-10 / (2 ** 0.5), -10 / (2 ** 0.5), -25),
    (10 / (2 ** 0.5), 10 / (2 ** 0.5), 25),
    only_fissionable=True,
)
sou.energy = openmc.stats.Discrete([1.0e9], [1.0])
sets = openmc.Settings()
sets.batches = 200
sets.inactive = 50
sets.particles = 16 * 625
sets.run_mode = "eigenvalue"
sets.keff_trigger = {"type": "std_dev", "threshold": 5e-4}
sets.trigger_active = True
sets.trigger_max_batches = 200
sets.sou = sou


# In[ ]:


tals = openmc.Tallies()
mesh = openmc.RegularMesh()
mesh.dimension = [100, 100]
mesh.lower_left = [-320, -320]
mesh.upper_right = [320, 320]
mesh_fil = openmc.MeshFilter(mesh)
tal = openmc.Tally(name="flux")
tal.filters = [mesh_fil]
tal.scores = ["flux", "fission"]
tals.append(tal)


# In[ ]:


mats.export_to_xml()
geo.export_to_xml()
sets.export_to_xml()
tals.export_to_xml()


# In[ ]:


openmc.run(output=False)


# In[ ]:


# sp = openmc.StatePoint("statepoint.200.h5")
# energy_bins = np.logspace(3, 9)
# probability, bin_edges = np.histogram(sp.source["E"], energy_bins, density=True)
# print(sum(probability * np.diff(energy_bins)))
# plt.semilogx(energy_bins[:-1], probability * np.diff(energy_bins), drawstyle="steps")
# plt.xlabel("Energy (eV)")
# plt.ylabel("Probability/eV")


# In[ ]:


# plt.quiver(
#     sp.source["r"]["x"],
#     sp.source["r"]["y"],
#     sp.source["u"]["x"],
#     sp.source["u"]["y"],
#     np.log(sp.source["E"]),
#     cmap="jet",
#     scale=20.0,
# )
# plt.colorbar()


# In[ ]:


# sp.keff


# In[ ]:


# sp.close()


# In[ ]:


# mod = openmc.Model(geometry=geo, settings=sets)
# operator = openmc.deplete.Operator(mod, path)
# timesteps = [140] * 16
# power = 1500 * 1e6
# timestep_units = "d"
# integrator = openmc.deplete.PredictorIntegrator(
#     operator=operator,
#     timesteps=timesteps,
#     power=power,
#     timestep_units=timestep_units,
# )
# integrator.integrate()


# In[ ]:


# results = openmc.deplete.Results("depletion_results.h5")
# time, k = results.get_keff()
# time /= 24 * 60 * 60
# # pyplot.errorbar(time, k[:, 0], yerr=k[:, 1])
# pyplot.plot(time, k[:, 0])
# pyplot.xlabel("Time [d]")
# pyplot.ylabel("$k_{eff}\pm \sigma$")


# In[ ]:


# import os
# import numpy as np

# for i in np.linspace(0.09, 0.12, 61, endpoint=True):
#     os.system("python homework_openmc.py {}".format(i))


# In[ ]:


import os 
sp = openmc.StatePoint('statepoint.200.h5')
k = sp.keff
sp.close()
try: 
    os.system('rm statepoint.*')
    os.system('rm particle*')
    os.system('echo "P(u):{} keff:{}" >> keff.txt'.format(x,k))
except Exception as e:
    print(e)


# In[ ]:


# import re
# import matplotlib.pyplot as plt

# f = open("keff.txt")
# f = f.read().split(sep="\n")
# x = []
# y = []
# for i in f[:-1]:
#     x.append(float(re.findall(pattern=r"P\(u\):(.*) keff", string=i)[0]))
#     y.append(float(re.findall(pattern=r"keff:(.*)\+", string=i)[0]))
# plt.figure(dpi=160)
# plt.plot(x, y)
# plt.grid()
# plt.xlabel("$P_{(U_{233})}$")
# plt.ylabel("$keff$")
# plt.title("$P_{(U_{233})}$-$keff$")

