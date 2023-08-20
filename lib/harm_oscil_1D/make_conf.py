import numpy as np
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
N_conf = 20000
N_dump = 5000
N_skip = 150
delta = 1


def weight(u, t, a):
    u_b_t = u[(t-1+u.shape[0]) % u.shape[0]]
    u_f_t = u[(t+1) % u.shape[0]]
    u_t = u[t]
    S = (u_f_t-u_t)**2/(2.0*a)+a*(u_f_t**2+u_t**2)/4.0 + \
        (u_t-u_b_t)**2/(2.0*a)+a*(u_t**2+u_b_t**2)/4.0
    return np.exp(-1.0*S)


def genr_conf(Un, delta, a):
    old_u = np.zeros(Un.shape[1])
    new_u = np.zeros(Un.shape[1])
    for n in range(Un.shape[0]):
        t_iter = np.array(range(Un.shape[1]))
        np.random.shuffle(t_iter)
        for t in t_iter:
            new_u[t] = old_u[t]+delta*np.random.uniform(-1, 1)
            if np.random.uniform(0, 1) < weight(new_u, t, a)/weight(old_u, t, a):
                old_u[t] = new_u[t]
            else:
                new_u[t] = old_u[t]
        Un[n, :] = new_u
    return Un


def save_conf(Un, N_dump, N_skip, a):
    Um = Un[N_dump:Un.shape[0]:N_skip, :]
    np.save("Um_"+str(N_dump)+"_"+str(N_skip)+"_" +
            str(Un.shape[0])+"_"+str(Un.shape[1])+"_"+str(a)+".npy", Um)
    return Um


def conf(N_conf, N_dump, N_skip, N_t, a, delta):
    Un = np.zeros((N_conf, N_t))
    Un = genr_conf(Un, delta, a)
    return save_conf(Un, N_dump, N_skip, a)


i = rank+16*0
a = 0.01+i*0.01
Um = conf(N_conf, N_dump, N_skip, 500, a, delta)
N_t = 20+i*20
Um = conf(N_conf, N_dump, N_skip, N_t, 0.01, delta)

i = rank+16*1
a = 0.01+i*0.01
Um = conf(N_conf, N_dump, N_skip, 500, a, delta)
N_t = 20+i*20
Um = conf(N_conf, N_dump, N_skip, N_t, 0.01, delta)

i = rank+16*2
a = 0.01+i*0.01
Um = conf(N_conf, N_dump, N_skip, 500, a, delta)
N_t = 20+i*20
Um = conf(N_conf, N_dump, N_skip, N_t, 0.01, delta)

i = rank+16*3
a = 0.01+i*0.01
Um = conf(N_conf, N_dump, N_skip, 500, a, delta)
N_t = 20+i*20
Um = conf(N_conf, N_dump, N_skip, N_t, 0.01, delta)
