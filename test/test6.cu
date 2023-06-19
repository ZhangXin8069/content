#include <iostream>
#include <random>
#include <cmath>
#include <chrono>
#include <cuda.h>

class Lattice
{
public:
    Lattice(int size) : L(size)
    {
        // Allocate memory for the lattice
        sites.resize(4 * L * L * L * L);
    }
    ~Lattice() {}
    // Accessor function to retrieve the value of a site
    __host__ double &operator()(int x, int y, int z, int t, int mu)
    {
        return sites[mu + 4 * (x + L * (y + L * (z + L * t)))];
    }
    // Accessor function to retrieve pointer to the lattice data
    __host__ double *data() { return sites.data(); }

private:
    int L;
    std::vector<double> sites;
};

__global__ void initializeLattice(double *lattice, int L)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < 4 * L * L * L * L)
    {
        lattice[id] = 1.0;
    }
}

__global__ void computeDslash(double *phi, double *gauge, double *dslash_phi, int L)
{
    __shared__ double shared_gauge[4 * 8 * 8];
    int mu = threadIdx.x / 64;
    int x = (blockIdx.x * 8 + (threadIdx.x / 8) % 8) % L;
    int y = (blockIdx.x * 8 + threadIdx.x / 64) % L;
    int z = (blockIdx.x / (L / 8) * 8 + threadIdx.x % 8) % L;
    int t = (blockIdx.x / (L / 8) / 8 * 8 + threadIdx.x / 512) % L;
    double tmp = 0.0;
    // Copy the gauge field to shared memory
    for (int nu = 0; nu < 4; nu++)
    {
        shared_gauge[nu + 4 * (threadIdx.x % 64 + threadIdx.x / 256 * 64)] = gauge[nu + 4 * (x + L * (y + L * (z + L * t)))];
        __syncthreads();
    }
    // Compute the forward and backward differences
    for (int nu = 0; nu < 4; nu++)
    {
        if (nu != mu)
        {
            int forward_id = nu + 4 * (x + L * (y + L * (z + L * t)));
            int backward_id = nu + 4 * (((x + L - 1) % L) + L * (((y + L - 1) % L) + L * (((z + L - 1) % L) + L * ((t + L - 1) % L))));
            tmp += shared_gauge[nu + 4 * (threadIdx.x % 64 + threadIdx.x / 256 * 64)] * phi[backward_id];
            tmp -= shared_gauge[nu + 4 * ((threadIdx.x + 64) % 256 + threadIdx.x / 256 * 64)] * phi[forward_id];
            __syncthreads();
        }
    }
    dslash_phi[mu + 4 * (x + L * (y + L * (z + L * t)))] = 2.0 * phi[mu + 4 * (x + L * (y + L * (z + L * t)))] - tmp;
}

__device__ void computeDaggerDslash(double *phi, double *gauge, double *dagger_dslash_phi, int L)
{
    __shared__ double shared_gauge[4 * 8 * 8];
    int mu = threadIdx.x / 64;
    int x = (blockIdx.x * 8 + (threadIdx.x / 8) % 8) % L;
    int y = (blockIdx.x * 8 + threadIdx.x / 64) % L;
    int z = (blockIdx.x / (L / 8) * 8 + threadIdx.x % 8) % L;
    int t = (blockIdx.x / (L / 8) / 8 * 8 + threadIdx.x / 512) % L;
    double tmp = 0.0;
    // Copy the gauge field to shared memory
    for (int nu = 0; nu < 4; nu++)
    {
        shared_gauge[nu + 4 * (threadIdx.x % 64 + threadIdx.x / 256 * 64)] = gauge[nu + 4 * (x + L * (y + L * (z + L * t)))];
        __syncthreads();
    }
    // Compute the forward and backward differences
    for (int nu = 0; nu < 4; nu++)
    {
        if (nu != mu)
        {
            int forward_id = nu + 4 * (x + L * (y + L * (z + L * t)));
            int backward_id = nu + 4 * (((x + L - 1) % L) + L * (((y + L - 1) % L) + L * (((z + L - 1) % L) + L * ((t + L - 1) % L))));
            tmp += shared_gauge[nu + 4 * (threadIdx.x % 64 + threadIdx.x / 256 * 64)] * phi[backward_id];
            tmp -= shared_gauge[nu + 4 * ((threadIdx.x + 64) % 256 + threadIdx.x / 256 * 64)] * phi[forward_id];
            __syncthreads();
        }
    }
    dagger_dslash_phi[mu + 4 * (x + L * (y + L * (z + L * t)))] = 2.0 * phi[mu + 4 * (x + L * (y + L * (z + L * t)))] - tmp;
}

__global__ void computeR(double *r, double *b, double *a, double *x, double *gauge, int L)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < 4 * L * L * L * L)
    {
        r[id] = b[id] - a[id] * x[id];
        computeDaggerDslash(r, gauge, r, L);
    }
}

__global__ void computeP(double *p, double *r, double *beta, double *p_ap, int L)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < 4 * L * L * L * L)
    {
        p[id] = r[id] + beta[0] * p_ap[id];
    }
}

__global__ void computeAlpha(double *alpha, double *r, double *p)
{
    __shared__ double shared_r[4 * 256];
    __shared__ double shared_p[4 * 256];
    int id = threadIdx.x;
    shared_r[id] = r[id];
    shared_p[id] = p[id];
    __syncthreads();
    double rho = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        rho += shared_r[mu + 4 * id] * shared_r[mu + 4 * id];
    }
    double pAp = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        pAp += shared_p[mu + 4 * id] * shared_p[mu + 4 * id];
    }
    alpha[0] = rho / pAp;
}

__global__ void computeAx(double *Ax, double *p, double *gauge, int L)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < 4 * L * L * L * L)
    {
        computeDaggerDslash(p, gauge, Ax, L);
    }
}

__global__ void computeBeta(double *beta, double *r, double *r_new)
{
    __shared__ double shared_r[4 * 256];
    __shared__ double shared_r_new[4 * 256];
    int id = threadIdx.x;
    shared_r[id] = r[id];
    shared_r_new[id] = r_new[id];
    __syncthreads();
    double rho_new = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        rho_new += shared_r_new[mu + 4 * id] * shared_r_new[mu + 4 * id];
    }
    double rho = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        rho += shared_r[mu + 4 * id] * shared_r[mu + 4 * id];
    }
    beta[0] = rho_new / rho;
}

__global__ void updateX(double *x, double *alpha, double *p, int L)
{
    int id = blockIdx.x * blockDim.x + threadIdx.x;
    if (id < 4 * L * L * L * L)
    {
        x[id] += alpha[0] * p[id];
    }
}

void cg(Lattice &x, Lattice &b, Lattice &gauge, int L, int max_it)
{
    int num_threads = 256;
    int num_blocks = (4 * L * L * L * L + num_threads - 1) / num_threads;
    double *d_x;
    double *d_b;
    double *d_gauge;
    double *d_r;
    double *d_p;
    double *d_Ap;
    double *d_beta;
    double *d_alpha;
    double *d_pAp;
    cudaMalloc(&d_x, 4 * L * L * L * L * sizeof(double));
    cudaMalloc(&d_b, 4 * L * L * L * L * sizeof(double));
    cudaMalloc(&d_gauge, 4 * L * L * L * L * sizeof(double));
    cudaMalloc(&d_r, 4 * L * L * L * L * sizeof(double));
    cudaMalloc(&d_p, 4 * L * L * L * L * sizeof(double));
    cudaMalloc(&d_Ap, 4 * L * L * L * L * sizeof(double));
    cudaMalloc(&d_beta, sizeof(double));
    cudaMalloc(&d_alpha, sizeof(double));
    cudaMalloc(&d_pAp, sizeof(double));
    cudaMemcpy(d_x, x.data(), 4 * L * L * L * L * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_b, b.data(), 4 * L * L * L * L * sizeof(double), cudaMemcpyHostToDevice);
    cudaMemcpy(d_gauge, gauge.data(), 4 * L * L * L * L * sizeof(double), cudaMemcpyHostToDevice);
    // Initialize the lattice to 1.0
    initializeLattice<<<num_blocks, num_threads>>>(d_x, L);
    // Initialize the residual
    computeR<<<num_blocks, num_threads>>>(d_r, d_b, d_Ap, d_x, d_gauge, L);
    cudaMemcpy(d_p, d_r, 4 * L * L * L * L * sizeof(double), cudaMemcpyDeviceToDevice);
    double pAp, beta;
    for (int i = 0; i < max_it; i++)
    {
        // Compute Ap
        computeAx<<<num_blocks, num_threads>>>(d_Ap, d_p, d_gauge, L);
        // Compute alpha
        computeAlpha<<<1, num_threads>>>(d_alpha, d_r, d_p);
        // Update x
        updateX<<<num_blocks, num_threads>>>(d_x, d_alpha, d_p, L);
        // Update r
        updateX<<<num_blocks, num_threads>>>(d_r, d_alpha, d_Ap, L);
        // Compute beta
        computeBeta<<<1, num_threads>>>(d_beta, d_r, d_p);
        // Update p
        computeP<<<num_blocks, num_threads>>>(d_p, d_r, d_beta, d_pAp, L);
        pAp = 0.0;
        // Compute pAp
        for (int mu = 0; mu < 4; mu++)
        {
            computeDslash<<<num_blocks, num_threads>>>(d_Ap + 4 * (L - 1) + mu * 4 * L * L * L, d_gauge, d_Ap + mu * 4 * L * L * L, L);
            pAp += d_p[mu * (L * L * L * L - 1)] * d_Ap[mu * (L * L * L * L - 1)];
        }
        cudaMemcpy(&beta, d_beta, sizeof(double), cudaMemcpyDeviceToHost);
        cudaMemcpy(&pAp, d_pAp, sizeof(double), cudaMemcpyDeviceToHost);
        std::cout << "Iteration " << i + 1 << ": ||r|| = " << sqrt(pAp) << std::endl;
        if (sqrt(pAp) < 1e-6)
            break;
    }
    cudaMemcpy(x.data(), d_x, 4 * L * L * L * L * sizeof(double), cudaMemcpyDeviceToHost);
    cudaFree(d_x);
    cudaFree(d_b);
    cudaFree(d_gauge);
    cudaFree(d_r);
    cudaFree(d_p);
    cudaFree(d_Ap);
    cudaFree(d_beta);
    cudaFree(d_alpha);
    cudaFree(d_pAp);
}

int main()
{
    int L = 8;
    int max_it = 1000;
    Lattice x(L);
    Lattice b(L);
    Lattice gauge(L);
    cg(x, b, gauge, L, max_it);
    return 0;
}