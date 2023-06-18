from re import I
import numpy as np
import taichi as ti
ti.init()


@ti.data_oriented
class biCGstab:
    def __init__(self, n: int, A: np.array, b: np.array, tol: float, limit: int):
        self.n = n
        self.A = A
        self.b = b
        self.tol = tol
        self.limit = limit

   
    import numpy as np
    def CG(self):
        x0 = np.random.rand(n, 1)
        x = x0
        converged = np.array([False])
        v = r = r0_hat = p = s = t = np.zeros((self.n, 1))
        r[:] = self.b[:]-np.dot(self.A, x)
        r0_hat[:] = r[:]
        rho0 = alpha = w = v[:] = p[:] = 1.0
        rho1 = np.dot(np.transpose(r0_hat), r)
        iters = 0
        while(True):
            print(iters,x)
            iters = iters+1
            converged[:] = (np.linalg.norm(r) < self.tol*np.linalg.norm(self.b))
            if converged == True or iters == self.limit:
                print(iters,":end")
                return x0, x, iters
            beta = (rho1/rho0)*(alpha/w)
            p[:, 0] = r[:, 0]+beta*(p[:, 0]-w*v[:, 0])
            v[:] = np.dot(self.A, p)
            alpha = rho1/np.dot(np.transpose(r0_hat), v)
            s[:, 0] = r[:, 0]-alpha*v[:, 0]
            t[:] = np.dot(self.A, s)
            w = np.dot(np.transpose(t), s)/np.dot(np.transpose(t), t)
            rho0 = rho1
            rho1 = -w*np.dot(np.transpose(r0_hat), t)
            # print(alpha, p, w, s)
            x[:, 0] = x[:, 0]+alpha*p[:, 0]+w*s[:, 0]
            r[:, 0] = s[:, 0]-w*t[:, 0]
            if iters<5:
                print(x[:,0])

    def report(self):
        x0, x, iters = self.CG()
        print('系数矩阵:A')
        print(self.A[:])
        print('右手边向量:b', self.b[:, 0])
        print('初始猜测值:x0', x0[:, 0])
        print('迭代次数', iters)
        print('解向量', x[:, 0])


if __name__ == "__main__":

    import numpy as np
    n = 1
    A = np.random.rand(n, n)
    b = np.random.rand(n, 1)
    # print(A,b)
    # A = np.array([[10, 1, -5], [-20, 3, 20], [5, 3, 5]])
    # b = np.array([[1], [2], [6]])
    # print(A,b)
    tol = 1e-6
    limit = 1e6
    begin = biCGstab(n=n, A=A, tol=tol, limit=limit, b=b)
    begin.report()
