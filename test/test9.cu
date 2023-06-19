#include <iostream>
#include <random>
#include <chrono>
#include <cmath>
#include <cuda_runtime.h>

// Forward declarations
class Lattice;
class FermiField;
class GaugeField;
class GaugeFieldGradient;
class GaugeFieldUpdate;
class StaggeredFermionOperator;
class ConjugateGradientSolver;

// Global constants
constexpr int N = 32; // Lattice size (NxNxNxN)
constexpr int V = N * N * N * N;
constexpr int GSIZE = 4 * V;
constexpr int FSIZE = 2 * V;
constexpr int BLOCK_DIM = 256;

// Global variables
__constant__ int dev_N;
__constant__ int dev_V;
__constant__ int dev_GSIZE;
__constant__ int dev_FSIZE;
// Host functions

// Initialize random number generator
std::default_random_engine generator(std::chrono::system_clock::now().time_since_epoch().count());
std::uniform_real_distribution<double> distribution(-0.5, 0.5);
auto rnd = std::bind(distribution, generator);

// Initialize gauge field with random numbers
void init_gauge_field(double *gauge)
{
    for (int i = 0; i < GSIZE; i++)
    {
        gauge[i] = rnd();
    }
}

// Initialize fermi field with random numbers
void init_fermi_field(double *fermi)
{
    for (int i = 0; i < FSIZE; i++)
    {
        fermi[i] = rnd();
    }
}

// Device functions

// Compute gauge field gradient kernel
__global__ void gauge_field_gradient_kernel(const double *__restrict__ gauge,
                                            double *__restrict__ gauge_grad)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dev_GSIZE)
    {
        const int t = idx / (dev_N * dev_N * dev_N);
        const int z = (idx % (dev_N * dev_N * dev_N)) / (dev_N * dev_N);
        const int y = (idx % (dev_N * dev_N)) / dev_N;
        const int x = idx % dev_N;
        const int idx_xp = ((t * dev_N + z) * dev_N + y) * dev_N + (x + 1);
        const int idx_xm = ((t * dev_N + z) * dev_N + y) * dev_N + (x - 1);
        const int idx_yp = ((t * dev_N + z) * dev_N + (y + 1)) * dev_N + x;
        const int idx_ym = ((t * dev_N + z) * dev_N + (y - 1)) * dev_N + x;
        const int idx_zp = ((t * dev_N + (z + 1)) * dev_N + y) * dev_N + x;
        const int idx_zm = ((t * dev_N + (z - 1)) * dev_N + y) * dev_N + x;
        gauge_grad[idx] = (gauge[idx_xp] - gauge[idx_xm] +
                           gauge[idx_yp] - gauge[idx_ym] +
                           gauge[idx_zp] - gauge[idx_zm]) /
                          2.0;
    }
}

// Update gauge field kernel
__global__ void gauge_field_update_kernel(const double *__restrict__ gauge,
                                          const double *__restrict__ gauge_grad,
                                          const double delta,
                                          double *__restrict__ new_gauge)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dev_GSIZE)
    {
        new_gauge[idx] = gauge[idx] - delta * gauge_grad[idx];
    }
}

// Staggered fermion operator kernel
__global__ void staggered_fermion_op_kernel(const double *__restrict__ fermi_in,
                                            const double *__restrict__ gauge,
                                            double *__restrict__ fermi_out)
{
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;

    if (idx < dev_FSIZE)
    {
        const int x = idx % dev_N;
        const int y = (idx % (dev_N * dev_N)) / dev_N;
        const int z = (idx % (dev_N * dev_N * dev_N)) / (dev_N * dev_N);
        const int t = idx / (dev_N * dev_N * dev_N);
        const int s = (x + y + z + t) % 2;
        const int sign = s ? -1 : 1;
        const int idx_xm = ((t * dev_N + z) * dev_N + (y - sign)) * dev_N + x; //
        const int idx_xp = ((t * dev_N + z) * dev_N + (y + sign)) * dev_N + x; //
        const int idx_ym = ((t * dev_N + (z - sign)) * dev_N + y) * dev_N + x;
        const int idx_yp = ((t * dev_N + (z + sign)) * dev_N + y) * dev_N + x;
        const int idx_zm = ((t * dev_N + z) * dev_N + (y - sign)) * dev_N + x;
        const int idx_zp = ((t * dev_N + z) * dev_N + (y + sign)) * dev_N + x;
        const int idx = idx_xp + (dev_N * dev_N * dev_N) * (s * dev_V) + (idx_xm - idx_xp) * (idx_xp >= s * dev_V);
        fermi_out[idx] = fermi_in[idx_xp] - sign * gauge[idx] * fermi_in[idx_ym];
        fermi_out[idx + dev_V] = fermi_in[idx_yp] + sign * gauge[idx] * fermi_in[idx_zm];
    }
}

// Conjugate gradient solver class
class ConjugateGradientSolver
{
public:
    ConjugateGradientSolver(const Lattice *lat, double epsilon, int max_iterations)
        : lat_(lat), epsilon_(epsilon), max_iterations_(max_iterations) {}

    // Solve linear system A x = b using the Conjugate Gradient method
    void solve(const double *b, double *x)
    {
        double *r = lat_->Tmp();
        double *p = lat_->Tmp() + dev_FSIZE;
        double *Ap = lat_->Tmp() + 2 * dev_FSIZE;
        double *tmp = lat_->Tmp() + 3 * dev_FSIZE;

        // Initialize variables
        double alpha = 0.0;
        double beta = 0.0;
        double rho = 0.0;
        double rho_old = 0.0;

        // Compute initial residual r = b - A x
        staggered_fermion_op(x, lat_->Gauge(), tmp);
        for (int i = 0; i < dev_FSIZE; i++)
        {
            r[i] = b[i] - tmp[i];
        }

        // Initialize p = r
        std::memcpy(p, r, dev_FSIZE * sizeof(double));

        // Compute initial rho = r^T r
        rho = dot_product(r, r);

        // Iterate until convergence or max_iterations reached
        for (int k = 0; k < max_iterations_; k++)
        {

            // Compute A p
            staggered_fermion_op(p, lat_->Gauge(), Ap);

            // Compute alpha = rho / p^T A p
            alpha = rho / dot_product(p, Ap);

            // Update solution x = x + alpha p
            for (int i = 0; i < dev_FSIZE; i++)
            {
                x[i] += alpha * p[i];
            }

            // Update residual r = r - alpha A p
            for (int i = 0; i < dev_FSIZE; i++)
            {
                r[i] -= alpha * Ap[i];
            }

            // Compute new rho = r^T r
            rho_old = rho;
            rho = dot_product(r, r);

            // Check for convergence
            if (sqrt(rho) < epsilon_)
            {
                break;
            }

            // Compute beta = rho / rho_old
            beta = rho / rho_old;

            // Update p = r + beta p
            for (int i = 0; i < dev_FSIZE; i++)
            {
                p[i] = r[i] + beta * p[i];
            }
        }
    }

private:
    const Lattice *lat_;
    double epsilon_;
    int max_iterations_;

    // Compute dot product of two vectors
    double dot_product(const double *x, const double *y)
    {
        double result = 0.0;
        for (int i = 0; i < dev_FSIZE; i++)
        {
            result += x[i] * y[i];
        }
        return result;
    }
};

// Lattice class
class Lattice
{
public:
    Lattice()
    {
        // Initialize device constants
        cudaMalloc(&dev_N, sizeof(int));
        cudaMalloc(&dev_V, sizeof(int));
        cudaMalloc(&dev_GSIZE, sizeof(int));
        cudaMalloc(&dev_FSIZE, sizeof(int));
        cudaMemcpy(dev_N, &N, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_V, &V, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_GSIZE, &GSIZE, sizeof(int), cudaMemcpyHostToDevice);
        cudaMemcpy(dev_FSIZE, &FSIZE, sizeof(int), cudaMemcpyHostToDevice);

        // Allocate device memory
        cudaMalloc(&gauge_, Gsizeof(double));
        cudaMalloc(&fermi_, FSIZE * sizeof(double));
        cudaMalloc(&tmp_, 4 * FSIZE * sizeof(double));

        // Initialize gauge and fermi fields
        init_gauge_field(gauge_);
        init_fermi_field(fermi_);
    }

    ~Lattice()
    {
        // Free device memory
        cudaFree(gauge_);
        cudaFree(fermi_);
        cudaFree(tmp_);
        cudaFree(dev_N);
        cudaFree(dev_V);
        cudaFree(dev_GSIZE);
        cudaFree(dev_FSIZE);
    }

    // Getter functions
    const double *Gauge() const { return gauge_; }
    const double *Fermi() const { return fermi_; }
    double *Tmp() const { return tmp_; }

    // Update gauge field using the Stable Conjugate Gradient algorithm
    void UpdateGauge(double delta)
    {
        GaugeFieldGradient gauge_grad(gauge_);
        GaugeFieldUpdate gauge_update(gauge_grad, delta);
        ConjugateGradientSolver cg_solver(this, 1e-10, 1000);
        cg_solver.solve(gauge_update.B(), gauge_update.X());
        gauge_update.Apply(gauge_);
    }

    // Apply staggered fermion operator
    void staggered_fermion_op(const double *fermi_in, const double *gauge, double *fermi_out)
    {
        StaggeredFermionOperator fermion_op(gauge);
        fermion_op.apply(fermi_in, fermi_out);
    }

private:
    double *gauge_;
    double *fermi_;
    double *tmp_;
};

// Gauge field gradient class
class GaugeFieldGradient
{
public:
    GaugeFieldGradient(const double *gauge)
    {
        cudaMalloc(&gauge_grad_, GSIZE * sizeof(double));
        gauge_field_gradient_kernel<<<(GSIZE + BLOCK_DIM - 1) / BLOCK_DIM, BLOCK_DIM>>>(gauge, gauge_grad_);
        cudaDeviceSynchronize();
    }

    ~GaugeFieldGradient()
    {
        cudaFree(gauge_grad_);
    }

    const double *Grad() const { return gauge_grad_; }

private:
    double *gauge_grad_;
};

// Gauge field update class
class GaugeFieldUpdate
{
public:
    GaugeFieldUpdate(const GaugeFieldGradient &gauge_grad, double delta)
        : delta_(delta)
    {
        cudaMalloc(&B_, GSIZE * sizeof(double));
        cudaMalloc(&X_, GSIZE * sizeof(double));
        cudaMemcpy(B_, gauge_grad.Grad(), GSIZE * sizeof(double), cudaMemcpyDeviceToDevice);
        scale<<<(GSIZE + BLOCK_DIM - 1) / BLOCK_DIM, BLOCK_DIM>>>(B_, -delta_);
        cudaDeviceSynchronize();
    }

    ~GaugeFieldUpdate()
    {
        cudaFree(B_);
        cudaFree(X_);
    }

    const double *B() const { return B_; }
    double *X() const { return X_; }

    void Apply(double *gauge) const
    {
        gauge_field_update_kernel<<<(GSIZE + BLOCK_DIM - 1) / BLOCK_DIM, BLOCK_DIM>>>(gauge, B_, delta_, X_);
        cudaDeviceSynchronize();
        std::memcpy(gauge, X_, GSIZE * sizeof(double));
    }

private:
    double *B_;
    double *X_;
    double delta_;
};

// Staggered fermion operator class
class StaggeredFermionOperator
{
public:
    StaggeredFermionOperator(const double *gauge)
    {
        cudaMalloc(&gauge_, GSIZE * sizeof(double));
        cudaMemcpy(gauge_, gauge, GSIZE * sizeof(double), cudaMemcpyDeviceToDevice);
    }

    ~StaggeredFermionOperator()
    {
        cudaFree(gauge_);
    }

    void apply(const double *fermi_in, double *fermi_out) const
    {
        staggered_fermion_op_kernel<<<(FSIZE + BLOCK_DIM - 1) / BLOCK_DIM, BLOCK_DIM>>>(fermi_in, gauge_, fermi_out);
        cudaDeviceSynchronize();
    }

private:
    double *gauge_;
};

// Main function
int main()
{
    Lattice lat;

    // Update gauge field
    lat.UpdateGauge(0.1);

    // Apply staggered fermion operator
    double *fermi_in = lat.Fermi();
    double *fermi_out = lat.Tmp();
    init_fermi_field(fermi_in);
    lat.staggered_fermion_op(fermi_in, lat.Gauge(), fermi_out);

    // Output result
    std::cout << "Fermi field:\n";
    for (int i = 0; i < FSIZE; i++)
    {
        std::cout << fermi_out[i] << " ";
        if ((i + 1) % N == 0)
        {
            std::cout << "\n";
        }
    }

    return 0;
}