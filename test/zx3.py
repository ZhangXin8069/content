# 稳定双共轭梯度法
import numpy as np
n = 33
converged = np.array([False])
v = np.zeros((n, 1))
r = np.zeros((n, 1))
r0_hat = np.zeros((n, 1))
p = np.zeros((n, 1))
s = np.zeros((n, 1))
t = np.zeros((n, 1))
A = np.random.rand(n, n)
b = np.random.rand(n, 1)
x = np.ones((n, 1))
tol = 1e-6
limit = 1e6
print('系数矩阵 A:\n', A)
print('右边向量 b:\n', b)
print('初始猜测值 x0:\n', x)
r[:] = b[:]-np.dot(A, x)
r0_hat[:] = r[:]
rho0 = alpha = w = v[:] = p[:] = 1.0
rho1 = np.dot(np.transpose(r0_hat), r)
iters = 0
while(True):
    iters = iters+1
    converged[:] = (np.linalg.norm(r) < tol*np.linalg.norm(b))
    if converged == True or iters == limit:
        break
    # if iters < 5:
    #     print("先验值x({}):\n{}".format(iters, x))
    beta = (rho1/rho0)*(alpha/w)
    p[:, 0] = r[:, 0]+beta*(p[:, 0]-w*v[:, 0])
    v[:] = np.dot(A, p)
    alpha = rho1/np.dot(np.transpose(r0_hat), v)
    s[:, 0] = r[:, 0]-alpha*v[:, 0]
    t[:] = np.dot(A, s)
    w = np.dot(np.transpose(t), s)/np.dot(np.transpose(t), t)
    rho0 = rho1
    rho1 = -w*np.dot(np.transpose(r0_hat), t)
    x[:, 0] = x[:, 0]+alpha*p[:, 0]+w*s[:, 0]
    r[:, 0] = s[:, 0]-w*t[:, 0]

print('迭代次数 iters:\n', iters)
print('解向量 x:\n', x)
import numpy as np
n = 4
converged = np.array([False])
v = np.zeros((n, 1))
r = np.zeros((n, 1))
r0_hat = np.zeros((n, 1))
p = np.zeros((n, 1))
s = np.zeros((n, 1))
t = np.zeros((n, 1))
A = np.random.rand(n, n)
b = np.random.rand(n, 1)
x = np.ones((n, 1))
tol = 1e-6
limit = 1e6
print('系数矩阵 A:\n', A)
print('右边向量 b:\n', b)
print('初始猜测值 x0:\n', x)
r[:] = b[:]-np.dot(A, x)
r0_hat[:] = r[:]
rho0 = alpha = w = v[:] = p[:] = 1.0
rho1 = np.dot(np.transpose(r0_hat), r)
iters = 0
while(True):
    iters = iters+1
    converged[:] = (np.linalg.norm(r) < tol*np.linalg.norm(b))
    if converged == True or iters == limit:
        break
    # if iters < 5:
    #     print("先验值x({}):\n{}".format(iters, x))
    beta = (rho1/rho0)*(alpha/w)
    p[:, 0] = r[:, 0]+beta*(p[:, 0]-w*v[:, 0])
    v[:] = np.dot(A, p)
    alpha = rho1/np.dot(np.transpose(r0_hat), v)
    s[:, 0] = r[:, 0]-alpha*v[:, 0]
    t[:] = np.dot(A, s)
    w = np.dot(np.transpose(t), s)/np.dot(np.transpose(t), t)
    rho0 = rho1
    rho1 = -w*np.dot(np.transpose(r0_hat), t)
    x[:, 0] = x[:, 0]+alpha*p[:, 0]+w*s[:, 0]
    r[:, 0] = s[:, 0]-w*t[:, 0]

print('迭代次数 iters:\n', iters)
print('解向量 x:\n', x)