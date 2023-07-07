import os
import sys

# test_dir = os.path.dirname(os.path.abspath(__file__))
# sys.path.insert(0, os.path.join(test_dir, ".."))

import numpy as np

from pyquda import init

os.environ["QUDA_RESOURCE_PATH"] = ".cache"
init()

Lx, Ly, Lz, Lt = 16, 16, 16, 32
Nd, Ns, Nc = 4, 4, 3
latt_size = [Lx, Ly, Lz, Lt]

def gamma(i: int):
    gamma1 = np.array([[0, 0, 0, 0+1j], [0, 0, 0+1j, 0], [0, 0-1j, 0, 0], [0-1j, 0, 0, 0]], dtype=np.complex128)
    gamma2 = np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]], dtype=np.complex128)
    gamma3 = np.array([[0, 0, 0+1j, 0], [0, 0, 0, 0-1j], [0-1j, 0, 0, 0], [0, 0+1j, 0, 0]], dtype=np.complex128)
    gamma4 = np.array([[0, 0, 1, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 1, 0, 0]], dtype=np.complex128)
    e = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype = np.complex128)
    if i == 0:
        return gamma1
    elif i == 1:
        return gamma2
    elif i == 2:
        return gamma3
    elif i == 3:
        return gamma4
    else:
        print('error i')
        return np.array([[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=np.complex128)

def MyDslash(U, a, b, kappa = 1):
    first_item = a # 第一项：delta函数
    e = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]], dtype=np.complex128)
    # print('a = ', a[0,0,0,0])
    b = np.zeros(a.shape, dtype=np.complex128)
    # b = b + a

    # print('U[0, 0, 0, 0, 0] = \n', U[0, 0, 0, 0, 0])
    # print('======================')
    # print('U[1, 0, 0, 0, 0] = \n', U[1, 0, 0, 0, 0])
    # print('======================')
    # print('U[2, 0, 0, 0, 0] = \n', U[2, 0, 0, 0, 0])
    # print('======================')
    # print('U[3, 0, 0, 0, 0] = \n', U[3, 0, 0, 0, 0])
    # print('======================')
    # print('U[0, 1, 0, 0, 0] = \n', U[0, 1, 0, 0, 0])
    # print('======================')

    # calculate every matrix on b
    for x in range(0, Lx):
        for y in range(0, Ly):
            for z in range(0, Lz):
                for t in range(0, Lt):
                    for i in range(0, 4):
                        if i == 0:
                            product1 = np.tensordot((e-gamma(i)), U[i, t, z, y, x], axes = 0)
                            second_item1 = np.dot(product1.swapaxes(1,2).reshape(12, 12), a[t, z, y, (x+1)%Lx].reshape(12, 1))
                            product2 = np.tensordot((e+gamma(i)), np.matrix(U[i, t, z, y, (x+Lx-1)%Lx], dtype = np.complex128).H, axes = 0)
                            second_item2 = np.dot(product2.swapaxes(1,2).reshape(12, 12), a[t, z, y, (x+Lx-1)%Lx].reshape(12, 1))
                            b[t, z, y, x] += kappa * (second_item1 + second_item2).reshape(4, 3)
                        elif i == 1:
                            product1 = np.tensordot((e-gamma(i)), U[i, t, z, y, x], axes = 0)
                            second_item1 = np.dot(product1.swapaxes(1,2).reshape(12, 12), a[t, z, (y+1)%Ly, x].reshape(12, 1))
                            product2 = np.tensordot((e+gamma(i)), np.matrix(U[i, t, z, (y+Ly-1)%Ly, x], dtype = np.complex128).H, axes = 0)
                            second_item2 = np.dot(product2.swapaxes(1,2).reshape(12, 12), a[t, z, (y+Ly-1)%Ly, x].reshape(12, 1))
                            b[t, z, y, x] += kappa * (second_item1 + second_item2).reshape(4, 3)
                        elif i == 2:
                            product1 = np.tensordot((e-gamma(i)), U[i, t, z, y, x], axes = 0)
                            second_item1 = np.dot(product1.swapaxes(1,2).reshape(12, 12), a[t, (z+1)%Lz, y, x].reshape(12, 1))
                            product2 = np.tensordot((e+gamma(i)), np.matrix(U[i, t, (z+Lz-1)%Lz, y, x], dtype = np.complex128).H, axes = 0)
                            second_item2 = np.dot(product2.swapaxes(1,2).reshape(12, 12), a[t, (z+Lz-1)%Lz, y, x].reshape(12, 1))
                            b[t, z, y, x] += kappa * (second_item1 + second_item2).reshape(4, 3)
                        else:
                            product1 = np.tensordot((e-gamma(i)), U[i, t, z, y, x], axes = 0)
                            second_item1 = np.dot(product1.swapaxes(1,2).reshape(12, 12), a[(t+1)%Lt, z, y, x].reshape(12, 1))
                            product2 = np.tensordot((e+gamma(i)), np.matrix(U[i, (t+Lt-1)%Lt, z, y, x], dtype = np.complex128).H, axes = 0)
                            second_item2 = np.dot(product2.swapaxes(1,2).reshape(12, 12), a[(t+Lt-1)%Lt, z, y, x].reshape(12, 1))
                            b[t, z, y, x] += kappa * (second_item1 + second_item2).reshape(4, 3)
    return b

def myApplyDslash(Mp, p, U_seed):
    import cupy as cp
    from pyquda import core, quda
    from pyquda.enum_quda import QudaParity
    from pyquda.field import LatticeFermion
    from pyquda.utils import gauge_utils
    # Set parameters in Dslash and use m=-3.5 to make kappa=1
    dslash = core.getDslash(latt_size, -3.5, 0, 0, anti_periodic_t=False)

    # Generate gauge and then load it
    U = gauge_utils.gaussGauge(latt_size, U_seed)
    #dslash.loadGauge(U)
    U1 = U.lexico()

    # # Dslash a = b
    # quda.dslashQuda(b.even_ptr, a.odd_ptr, dslash.invert_param, QudaParity.QUDA_EVEN_PARITY)
    # quda.dslashQuda(b.odd_ptr, a.even_ptr, dslash.invert_param, QudaParity.QUDA_ODD_PARITY)

    # Save b to Mp
    # Mp[:] = b.lexico()

    Mp[:] = MyDslash(U1, p, Mp)
    # Return gauge as a ndarray with shape (Nd, Lt, Lz, Ly, Lx, Ns, Ns)
    return U1








def applyDslash(Mp, p, U_seed):
    import cupy as cp
    from pyquda import core, quda
    from pyquda.enum_quda import QudaParity
    from pyquda.field import LatticeFermion
    from pyquda.utils import gauge_utils

    # Set parameters in Dslash and use m=-3.5 to make kappa=1
    dslash = core.getDslash(latt_size, -3.5, 0, 0, anti_periodic_t=False)

    # Generate gauge and then load it
    U = gauge_utils.gaussGauge(latt_size, U_seed)
    dslash.loadGauge(U)

    # Load a from p and allocate b
    a = LatticeFermion(latt_size, cp.asarray(core.cb2(p, [0, 1, 2, 3])))
    b = LatticeFermion(latt_size)

    # Dslash a = b
    quda.dslashQuda(b.even_ptr, a.odd_ptr, dslash.invert_param, QudaParity.QUDA_EVEN_PARITY)
    quda.dslashQuda(b.odd_ptr, a.even_ptr, dslash.invert_param, QudaParity.QUDA_ODD_PARITY)

    # Save b to Mp
    Mp[:] = b.lexico()

    # Return gauge as a ndarray with shape (Nd, Lt, Lz, Ly, Lx, Ns, Ns)
    return U.lexico()


p = np.zeros((Lt, Lz, Ly, Lx, Ns, Nc), np.complex128)
p[0, 0, 0, 0, 0, 0] = 1
Mp = np.zeros((Lt, Lz, Ly, Lx, Ns, Nc), np.complex128)

U = applyDslash(Mp, p, 0)
print(Mp[0, 0, 0, 1])


# p1 = np.zeros((Lt, Lz, Ly, Lx, Ns, Nc), np.complex128)
# p[0, 0, 0, 0, 0, 0] = 1
Mp1 = np.zeros((Lt, Lz, Ly, Lx, Ns, Nc), np.complex128)

U = myApplyDslash(Mp1, p, 0)
print(Mp1[0, 0, 0, 1])
print(np.linalg.norm(Mp-Mp1))