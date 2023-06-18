#include <complex>
#include <iostream>
#include <mpi.h>
#include <random>
#include <vector>

template <typename T>
class BiCGStabSolver
{
public:
    // Constructor
    BiCGStabSolver(const std::vector<std::vector<T>> &A,
                   const std::vector<T> &b,
                   int max_iter = 1000,
                   double tolerance = 1e-6)
        : A_(A),
          b_(b),
          max_iter_(max_iter),
          tolerance_(tolerance),
          comm_(MPI_COMM_WORLD),
          rank_(0),
          size_(1)
    {
        // Initialize MPI
        MPI_Init(nullptr, nullptr);
        MPI_Comm_rank(comm_, &rank_);
        MPI_Comm_size(comm_, &size_);
    }

    // Destructor
    ~BiCGStabSolver() { MPI_Finalize(); }

    // Solve linear system of equations using the BiCGStab method
    std::vector<T> Solve()
    {
        // Get size of system
        int n = b_.size();

        // Initialize solution vector and residual
        std::vector<T> x(n, T{});
        std::vector<T> r = b_;

        // Initialize temporary vectors
        std::vector<T> r_hat(n);
        std::vector<T> p(n);
        std::vector<T> v(n);
        std::vector<T> s(n);
        std::vector<T> t(n);

        // Initialize variables
        T rho_prev = T{};
        T alpha = T{};
        T omega = T{};

        // Compute initial residual
        T rho = DotProduct(r, r);
        T rho_tol = tolerance_ * rho;

        // Initialize parity variables
        int parity = 0;
        int block_size = n / size_;
        int start = rank_ * block_size;
        int end = start + block_size;
        if (rank_ == size_ - 1)
        {
            end = n;
        }

        // Iterate until convergence or maximum iterations reached
        for (int i = 0; i < max_iter_; i++)
        {
            // Compute conjugate residual
            r_hat = Conjugate(r);

            // Compute parity-based BiCGStab method
            if (parity == 0)
            {
                // Compute rho
                rho = DotProduct(r, r_hat);

                // Check for convergence
                if (rho < rho_tol)
                {
                    break;
                }

                // Compute beta
                if (i == 0)
                {
                    p = r;
                }
                else
                {
                    T beta = (rho / rho_prev) * (alpha / omega);
                    p = r + beta * (p - omega);
                }

                // Compute v
                v = MatrixVectorProduct(A_, p, start, end);
                alpha = rho / DotProduct(r_hat, v);
                s = r - alpha * v;
            }
            else
            {
                // Compute t
                t = MatrixVectorProduct(A_, s, start, end);
                omega = DotProduct(t, s) / DotProduct(t, t);
                x += omega * s;
                r = s - omega * t;
            }

            // Update parity and rho_prev
            parity = (parity + 1) % 2;
            rho_prev = rho;
        }

        // Return solution vector
        return x;
    }

private:
    // Matrix-vector product
    std::vector<T> MatrixVectorProduct(const std::vector<std::vector<T>> &A,
                                       const std::vector<T> &x,
                                       int start,
                                       int end)
    {
        int n = A.size();
        std::vector<T> y(n, T{});
        for (int i = start; i < end; i++)
        {
            for (int j = 0; j < n; j++)
            {
                y[i] += A[i][j] * x[j];
            }
        }
        return y;
    }
    // Vector element-wise product
    std::vector<T> VectorElementwiseProduct(const std::vector<T> &x,
                                            const std::vector<T> &y)
    {
        int n = x.size();
        std::vector<T> z(n);
        for (int i = 0; i < n; i++)
        {
            z[i] = x[i] * y[i];
        }
        return z;
    }
    // Dot product
    T DotProduct(const std::vector<T> &x, const std::vector<T> &y)
    {
        int n = x.size();
        T sum = T{};
        for (int i = 0; i < n; i++)
        {
            sum += x[i] * y[i];
        }
        return sum;
    }

    // Conjugate vector
    std::vector<T> Conjugate(const std::vector<T> &x)
    {
        int n = x.size();
        std::vector<T> y(n);
        for (int i = 0; i < n; i++)
        {
            y[i] = std::conj(x[i]);
        }
        return y;
    }

    // Matrix
    std::vector<std::vector<T>> A_;

    // Right-hand side vector
    std::vector<T> b_;

    // Maximum number of iterations
    int max_iter_;

    // Tolerance for convergence
    double tolerance_;

    // MPI communicator
    MPI_Comm comm_;

    // MPI rank
    int rank_;

    // MPI size
    int size_;
};