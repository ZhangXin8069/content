#include <mpi.h>
#include <iostream>
#include <vector>

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
     int n = 8;

     // Create a matrix and a vector of the appropriate size.
     std::vector<std::vector<double>> A(n, std::vector<double>(n));
     std::vector<double> b(n);
     std::vector<double> x(n);

     // Initialize the matrix and the vector.
     for (int i = 0; i < n; i++)
     {
          for (int j = 0; j < n; j++)
          {
               A[i][j] = (i == j - 1 || i == j) ? 2 : 0;
          }
          b[i] = i + 1;
     }

     // Compute the local part of the solution vector.
     std::vector<double> local_x(n);
     for (int i = rank * 2; i < rank * 2 + 2; i++)
     {
          local_x[i] = b[i];
          for (int j = 0; j < n; j++)
          {
               if (j != i)
               {
                    local_x[i] -= A[i][j] * x[j];
               }
          }
          local_x[i] /= A[i][i];
     }

     // Asynchronously send and receive the local parts of the solution vector.
     MPI_Request send_request, recv_request;
     if (rank != 0)
     {
          MPI_Isend(&local_x[rank * 2], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &send_request);
          MPI_Irecv(&x[rank * 2 - 1], 1, MPI_DOUBLE, rank - 1, 0, MPI_COMM_WORLD, &recv_request);
     }
     if (rank != size - 1)
     {
          MPI_Isend(&local_x[rank * 2 + 1], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &send_request);
          MPI_Irecv(&x[rank * 2 + 2], 1, MPI_DOUBLE, rank + 1, 0, MPI_COMM_WORLD, &recv_request);
     }

     // Copy the local part of the solution vector into the global solution vector.

     for (int i = rank * 2; i < rank * 2 + 2; i++)
     {
          x[i] = local_x[i];
     }

     // Print the solution vector.
     if (rank == 0)
     {
          for (int i = 0; i < n; i++)
          {
               std::cout << x[i] << " ";
          }
          std::cout << std::endl;
     }

     // Finalize the MPI environment.
     MPI_Finalize();
}
