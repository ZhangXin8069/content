#include <iostream>
#include "my_vector.h"
int main()
{
    MPI_Init(NULL, NULL);
    int lattice_x(3 * 4 * 500), lattice_t(3 * 4 * 5), lattice_spin(2); // const lattice_spin=2
    int size(lattice_x * lattice_t * lattice_spin);
    int MAX_ITER(1e3);
    double TOL(1e-5), start, end;
    ComplexVector b(size);
    ComplexVector U(size);
    b.clean_1();
    b[0] = 10;
    U.clean_1();
    b.clean_rand();
    U.clean_rand();
    std::complex<double> mass(1);
    mass = b.dotX(U);
    start = MPI_Wtime();
    end = MPI_Wtime();
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    std::cout
        << "########"
        << "time cost:"
        << end - start
        << "s"
        << mass / double(size)
        << std::endl;
    MPI_Finalize();

    return 0;
}