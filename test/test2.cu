#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define L 4
#define N (L * L * L * L)
#define MAX_ITER 1000
#define TOL 1e-12

// Define lattice site indices
#define ix(i, j, k, t) ((t)*L * L * L + (k)*L * L + (j)*L + (i))

// Define lattice structure
class Lattice
{
public:
    double *fermion_field;
    double *gauge_field;

    // Constructor
    Lattice()
    {
        cudaMallocManaged(&fermion_field, N * sizeof(double));
        cudaMallocManaged(&gauge_field, 4 * N * sizeof(double));
        for (int i = 0; i < N; i++)
        {
            fermion_field[i] = 0.0;
            for (int mu = 0; mu < 4; mu++)
            {
                gauge_field[4 * i + mu] = 1.0;
            }
        }
    }

    // Destructor
    ~Lattice()
    {
        cudaFree(fermion_field);
        cudaFree(gauge_field);
    }
};

// Compute Wilson fermion operator on GPU
__device__ double wilson_fermion_op(Lattice *lat, int i, int j)
{
    int t0 = i / (L * L * L);
    int x0 = (i / L) % L;
    int y0 = (i % L) % L;
    int z0 = i % L;
    int t1 = j / (L * L * L);
    int x1 = (j / L) % L;
    int y1 = (j % L) % L;
    int z1 = j % L;
    double gamma[5][4] = {
        {0.0, 0.0, 0.0, 1.0},
        {0.0, 0.0, 1.0, 0.0},
        {0.0, 0.0, 0.0, -1.0},
        {0.0, -1.0, 0.0, 0.0},
        {1.0, 0.0, 0.0, 0.0}};
    double wf = 0.0;
    for (int mu = 0; mu < 4; mu++)
    {
        int i_mu = 4 * i + mu;
        int j_mu = 4 * j + mu;
        if (i_mu == j_mu)
        {
            continue;
        }
        int i_mu_p = 4 * ((i + (1 << mu)) % N) + mu;
        int i_mu_m = 4 * ((i - (1 << mu) + N) % N) + mu;
        double U_mu_p = lat->gauge_field[i_mu];
        double U_mu_m = lat->gauge_field[i_mu_m];
        double U_mu_p_inv = lat->gauge_field[i_mu_p];
        U_mu_p_inv = 1.0 / U_mu_p_inv;
        double phi = gamma[0][mu];
        phi *= U_mu_m * lat->fermion_field[j];
        phi -= U_mu_p_inv * lat->fermion_field[j_mu];
        wf += phi;
    }
    wf += 2.0 * lat->fermion_field[j];
    return wf;
}

// Compute gauge field gradient on GPU
__global__ void gauge_field_grad_kernel(Lattice *lat, double *grad)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= 4 * N)
    {
        return;
    }
    grad[i] = 0.0;
    int t = i / (4 * L * L * L);
    int x = (i / (4 * L * L)) % L;
    int y = (i / (4 * L)) % L;
    int z = (i / 4) % L;
    int mu = i % 4;
    for (int i0 = 0; i0 < L; i0++)
    {
        for (int i1 = 0; i1 < L; i1++)
        {
            for (int i2 = 0; i2 < L; i2++)
            {
                int j = ix(i0, i1, i2, t);
                int j_mu = 4 * j + mu;
                int j_mu_p = 4 * ((ix(i0, i1, i2, (t + 1) % L) * 4 + mu) % (4 * N)) + mu;
                int j_mu_m = 4 * ((ix(i0, i1, i2, (t - 1 + L) % L) * 4 + mu) % (4 * N)) + mu;
                int i_mu = 4 * i + mu;
                double U_mu = lat->gauge_field[i_mu];
                double U_mu_p = lat->gauge_field[j_mu_p];
                double U_mu_m = lat->gauge_field[j_mu_m];
                double U_mu_p_inv = 1.0 / U_mu_p;
                double U_mu_m_inv = 1.0 / U_mu_m;
                double phi = 0.0;
                phi += U_mu * wilson_fermion_op(lat, i, j_mu);
                phi -= U_mu_p_inv * wilson_fermion_op(lat, i, j_mu_p);
                phi -= U_mu_m_inv * wilson_fermion_op(lat, i, j_mu_m);
                grad[i_mu] += phi;
            }
        }
    }
}

// Compute gauge field gradient on CPU
void gauge_field_grad(Lattice *lat, double *grad)
{
    int block_size = 256;
    int num_blocks = (4 * N + block_size - 1) / block_size;
    gauge_field_grad_kernel<<<num_blocks, block_size>>>(lat, grad);
    cudaDeviceSynchronize();
}

// Initialize lattice
void init_lattice(Lattice *lat)
{
    lat->fermion_field = NULL;
    lat->gauge_field = NULL;
    cudaMallocManaged(&lat->fermion_field, N * sizeof(double));
    cudaMallocManaged(&lat->gauge_field, 4 * N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        lat->fermion_field[i] = 0.0;
        for (int mu = 0; mu < 4; mu++)
        {
            lat->gauge_field[4 * i + mu] = 1.0;
        }
    }
}

// Free lattice memory
void free_lattice(Lattice *lat)
{
    cudaFree(lat->fermion_field);
    cudaFree(lat->gauge_field);
}

// Compute LQCD stable conjugate gradient on GPU
void lqcd_cg(Lattice *lat, double *b, double *x)
{
    double *r = NULL;
    double *p = NULL;
    double *Ap = NULL;
    cudaMallocManaged(&r, N * sizeof(double));
    cudaMallocManaged(&p, N * sizeof(double));
    cudaMallocManaged(&Ap, N * sizeof(double));
    for (int i = 0; i < N; i++)
    {
        r[i] = b[i];
        p[i] = b[i];
    }
    double r_norm_sq = 0.0;
    double r_norm_sq_old = 0.0;
    double alpha = 0.0;
    double beta = 0.0;
    gauge_field_grad(lat, Ap);
    for (int i = 0; i < N; i++)
    {
        Ap[i] -= b[i];
        p[i] = r[i] - Ap[i];
        r_norm_sq += r[i] * r[i];
    }
    int iter = 0;
    while (iter < MAX_ITER)
    {
        gauge_field_grad(lat, Ap);
        double Ap_dot_p = 0.0;
        for (int i = 0; i < N; i++)
        {
            Ap_dot_p += Ap[i] * p[i];
        }
        alpha = r_norm_sq / Ap_dot_p;
        for (int i = 0; i < N; i++)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }
        r_norm_sq_old = r_norm_sq;
        r_norm_sq = 0.0;
        for (int i = 0; i < N; i++)
        {
            r_norm_sq += r[i] * r[i];
        }
        if (sqrt(r_norm_sq) < TOL)
        {
            break;
        }
        beta = r_norm_sq / r_norm_sq_old;
        for (int i = 0; i < N; i++)
        {
            p[i] = r[i] + beta * p[i];
        }
        iter++;
    }
    cudaFree(r);
    cudaFree(p);
    cudaFree(Ap);
}

int main()
{
    Lattice *lat = new Lattice();
    double *b = NULL;
    double *x = NULL;
    cudaMallocManaged(&b, N * sizeof(double));
    cudaMallocManaged(&x, N * sizeof(double));
    // Initialize source vector and solution vector
    for (int i = 0; i < N; i++)
    {
        b[i] = 1.0;
        x[i] = 0.0;
    }
    // Compute LQCD stable conjugate gradient
    lqcd_cg(lat, b, x);

    // Print solution vector
    for (int i = 0; i < N; i++)
    {
        printf("%f ", x[i]);
    }
    printf("\n");

    // Free memory
    cudaFree(b);
    cudaFree(x);
    delete lat;

    return 0;
}