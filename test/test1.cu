#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

// Lattice size
const int N = 16;

// Precision type
typedef float real_t;

// Lattice class
class Lattice
{
public:
    Lattice()
    {
        cudaMalloc(&fermi, N * N * N * N * sizeof(real_t));
        cudaMalloc(&gauge, N * N * N * N * sizeof(real_t));
        cudaMalloc(&gauge_grad, N * N * N * N * sizeof(real_t));
        cudaMalloc(&tmp, N * N * N * N * sizeof(real_t));
        cudaMemset(fermi, 0, N * N * N * N * sizeof(real_t));
        cudaMemset(gauge, 0, N * N * N * N * sizeof(real_t));
        cudaMemset(gauge_grad, 0, N * N * N * N * sizeof(real_t));
        cudaMemset(tmp, 0, N * N * N * N * sizeof(real_t));
    }

    ~Lattice()
    {
        cudaFree(fermi);
        cudaFree(gauge);
        cudaFree(gauge_grad);
        cudaFree(tmp);
    }

    real_t *Fermi() const { return fermi; }
    real_t *Gauge() const { return gauge; }
    real_t *GaugeGrad() const { return gauge_grad; }
    real_t *Tmp() const { return tmp; }

private:
    real_t *fermi;
    real_t *gauge;
    real_t *gauge_grad;
    real_t *tmp;
};

// Wilson fermion operator
__device__ void wilson_fermion_op(const real_t *__restrict__ fermi,
                                  const real_t *__restrict__ gauge,
                                  const int i, const int j, const int k, const int t,
                                  const Lattice *__restrict__ lat,
                                  real_t *__restrict__ tmp)
{
    const int idx = t * N * N * N + k * N * N + j * N + i;

    // Compute Wilson fermion operator
    tmp[idx] = fermi[idx] - 0.5 * (gauge[idx] * fermi[idx + 1] +
                                   gauge[idx - 1] * fermi[idx - 1] +
                                   gauge[idx + N] * fermi[idx + N] +
                                   gauge[idx - N] * fermi[idx - N] +
                                   gauge[idx + N * N] * fermi[idx + N * N] +
                                   gauge[idx - N * N] * fermi[idx - N * N] +
                                   gauge[idx + N * N * N] * fermi[idx + N * N * N] +
                                   gauge[idx - N * N * N] * fermi[idx - N * N * N]);
}

// Gauge field gradient kernel
__global__ void gauge_field_grad_kernel(const real_t *__restrict__ fermi,
                                        const real_t *__restrict__ gauge,
                                        const Lattice *__restrict__ lat,
                                        real_t *__restrict__ gauge_grad)
{
    const int tx = threadIdx.x;
    const int ty = threadIdx.y;
    const int tz = threadIdx.z;
    const int bx = blockIdx.x;
    const int by = blockIdx.y;
    const int bz = blockIdx.z;

    const int ix = tx + bx * blockDim.x;
    const int iy = ty + by * blockDim.y;
    const int iz = tz + bz * blockDim.z;

    if (ix >= N || iy >= N || iz >= N)
    {
        return;
    }

    for (int t = 0; t < N; t++)
    {
        const int idx = t * N * N * N + iz * N * N + iy * N + ix;

        // Compute gauge field gradient
        const int idx_xm = t * N * N * N + iz * N * N + iy * N + (ix - 1 + N) % N;
        const int idx_xp = t * N * N * N + iz * N * N + iy * N + (ix + 1) % N;
        const int idx_ym = t * N * N * N + iz * N * N + ((iy - 1 + N) % N) * N + ix;
        const int idx_yp = t * N * N * N + iz * N * N + ((iy + 1) % N) * N + ix;
        const int idx_zm = t * N * N * N + ((iz - 1 + N) % N) *N *N + iy *N + ix const int idx_zp = t * N * N * N + ((iz + 1) % N) * N * N + iy * N + ix;

        gauge_grad[idx] = 0.5 * (gauge[idx_xp] - gauge[idx_xm] + gauge[idx_yp] - gauge[idx_ym] + gauge[idx_zp] - gauge[idx_zm]) / lat->Tmp()[idx];
    }
}

// LQCD stable conjugate gradient algorithm
void lqcd_cg(const real_t *__restrict__ fermi,
             const real_t *__restrict__ gauge,
             const Lattice *__restrict__ lat,
             const int max_iters, const real_t tol)
{
    real_t *x = lat->Tmp();
    real_t *r = lat->GaugeGrad();
    real_t *p = lat->Fermi();
    real_t *Ap = lat->Gauge();

    // Compute initial residue
    cudaMemcpy(r, fermi, N * N * N * N * sizeof(real_t), cudaMemcpyDeviceToDevice);
    wilson_fermion_op(p, gauge, 0, 0, 0, 0, lat, Ap);
    for (int i = 0; i < N * N * N * N; i++)
    {
        r[i] -= Ap[i];
    }

    real_t rnorm1 = 0;
    for (int i = 0; i < N * N * N * N; i++)
    {
        rnorm1 += r[i] * r[i];
    }

    // Main iteration loop
    int iter = 0;
    while (iter < max_iters)
    {
        // Compute Ap
        wilson_fermion_op(p, gauge, 0, 0, 0, 0, lat, Ap);

        // Compute alpha
        real_t alpha = 0;
        for (int i = 0; i < N * N * N * N; i++)
        {
            alpha += p[i] * Ap[i];
        }
        alpha = rnorm1 / alpha;

        // Update x and r
        for (int i = 0; i < N * N * N * N; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        // Compute new residue norm
        real_t rnorm2 = 0;
        for (int i = 0; i < N * N * N * N; i++)
        {
            rnorm2 += r[i] * r[i];
        }

        // Check convergence
        if (rnorm2 < tol * tol)
        {
            break;
        }

        // Compute beta
        real_t beta = rnorm2 / rnorm1;

        // Update p
        for (int i = 0; i < N * N * N * N; i++)
        {
            p[i] = r[i] + beta * p[i];
        }

        rnorm1 = rnorm2;
        iter++;
    }

    printf("Converged in %d iterations.\n", iter);
}

int main()
{
    // Initialize lattice
    Lattice *lat = new Lattice();

    // Enter random fermi and gauge fields
    for (int i = 0; i < N * N * N * N; i++)
    {
        lat->Fermi()[i] = rand() / (real_t)RAND_MAX;
        lat->Gauge()[i] = rand() / (real_t)RAND_MAX;
    }

    // Compute gauge field gradient
    dim3 block(8, 8, 8);
    dim3 grid((N + block.x - 1) / block.x,
              (N + block.y - 1) / block.y,
              (N + block.z - 1) / block.z);
    gauge_field_grad_kernel<<<grid, block>>>(lat->Fermi(), lat->Gauge(), lat, lat->Tmp());

    // Run LQCD stable conjugate gradient algorithm
    const int max_iters = 1000;
    const real_t tol = 1e-6;
    lqcd_cg(lat->Fermi(), lat->Gauge(), lat, max_iters, tol);

    // Output results
    cudaMemcpy(lat->Tmp(), lat->Fermi(), N * N * N * N * sizeof(real_t), cudaMemcpyDeviceToDevice);
    wilson_fermion_op(lat->Tmp(), lat->Gauge(), 0, 0, 0, 0, lat, lat->Tmp());
    printf("Result: %f\n", lat->Tmp()[0]);

    // Clean up
    delete lat;

    return 0;
}