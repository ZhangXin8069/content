#include <complex>
#include <mpi.h>
#include <iostream>
#include <vector>
#include <random>

int main(int argc, char **argv)
{
    // Initialize the MPI environment.
    MPI_Init(NULL, NULL);

    // Get the number of processes.
    int size;
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Get the rank of the process.
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // Get the size of the matrix.
    int n = 2;

    // Create a matrix and a vector of the appropriate size.
    std::vector<std::vector<std::complex<double>>> A(n, std::vector<std::complex<double>>(n));
    std::vector<std::complex<double>> b(n);
    std::vector<std::complex<double>> x(n);

    // Initialize the random number generator.
    std::mt19937 generator(rank);
    std::uniform_real_distribution<double> distribution(-1.0, 1.0);

    // Randomly initialize the matrix and the vector.
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            A[i][j] = std::complex<double>(distribution(generator), distribution(generator));
        }
        b[i] = std::complex<double>(distribution(generator), distribution(generator));
    }

    // Print the matrix and the vector.
    std::cout << "Matrix:" << std::endl;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
        {
            std::cout << A[i][j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << "Vector:" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << b[i] << " ";
    }
    std::cout << std::endl;

    // Perform odd-even preconditioning on the matrix and the vector.
    if (rank % 2 == 0)
    {
        for (int i = rank * 2; i < rank * 2 + 2; i++)
        {
            if (i % 2 == 0)
            {
                b[i] /= A[i][i];
                for (int j = i + 1; j < n; j++)
                {
                    b[j] -= A[j][i] * b[i];
                }
            }
        }
    }
    else
    {
        for (int i = rank * 2; i < rank * 2 + 2; i++)
        {
            if (i % 2 == 1)
            {
                b[i] /= A[i][i];
                for (int j = 0; j < n; j++)
                {
                    A[i][j] /= A[i][i];
                }
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Initialize the biCGstab solver.
    std::vector<std::complex<double>> r(b);
    std::vector<std::complex<double>> p(b);

    // Enter the iteration loop.
    for (int k = 0; k < n; k++)
    {
        // Compute the value of alpha.
        std::complex<double> alpha = 0.0;
        for (int i = 0; i < n; i++)
        {
            alpha += std::conj(r[i]) * r[i];
        }

        // Compute the values of s and omega.
        std::vector<std::complex<double>> s(n);
        std::vector<std::complex<double>> omega(n);
        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                s[i] += A[i][j] * p[j];
                omega[i] += std::conj(A[j][i]) * s[j];
            }
        }
        std::complex<double> sigma = 0.0;
        for (int i = 0; i < n; i++)
        {
            sigma += std::conj(s[i]) * p[i];
        }
        std::complex<double> omega_alpha = 0.0;
        for (int i = 0; i < n; i++)
        {
            omega_alpha += std::conj(omega[i]) * r[i];
        }
        std::complex<double> beta = (alpha * sigma) / (omega_alpha * sigma);
        std::complex<double> gamma = (alpha * omega_alpha) / (omega_alpha * sigma);

        // Update the vectors x, r, and p.
        for (int i = 0; i < n; i++)
        {
            x[i] += beta * p[i] + gamma * s[i];
            r[i] -= beta * omega[i] + gamma * s[i];
            p[i] = r[i] + (gamma / beta) * p[i];
        }
    }

    // Print the solution.
    std::cout << "Solution:" << std::endl;
    for (int i = 0; i < n; i++)
    {
        std::cout << x[i].real() << " + " << x[i].imag() << "i ";
    }
    std::cout << std::endl;

    // Finalize the MPI environment.
    MPI_Finalize();
}
