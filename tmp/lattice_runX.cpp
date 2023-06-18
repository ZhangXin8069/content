#include "lattice_vec.hpp"
#include "lattice_dslashX.hpp"
#include "lattice_bicgX.hpp"
#include <iostream>
#include <ctime>
int main()
{
    MPI_Init(NULL, NULL);
    int lat_x(128);
    int lat_t(16*3*5);
    int lat_spin(2); // const lat_spin=2
    int size(lat_x * lat_t * lat_spin);
    int MAX_ITER(1e2);
    double TOL(1e-5);
    lattice_vec b(size);
    b.clean1();
    b[0] = 10.0;
    lattice_vec U(size);
    U.rand();
    double mass(1);
    bool dag(true);
    clock_t start = clock();
    lattice_dslash Dslash(U, lat_x, lat_t, lat_spin, mass, dag);
    lattice_bicg Bicg(MAX_ITER, TOL, b, Dslash);
    Bicg.solve();
    clock_t end = clock();
    std::cout
        << "################"
        << "time cost:"
        << (double)(end - start) / CLOCKS_PER_SEC
        << "s"
        << std::endl;
    // b.delete_();
    // U.delete_();
    MPI_Finalize();
    return 0;
}
