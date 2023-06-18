#include "lattice_vec.hpp"
#include "lattice_dslash.hpp"
#include "lattice_bicg.hpp"
#include <iostream>
#include <ctime>
int main()
{
    int lat_x(100);
    int lat_t(100);
    int lat_spin(2); // const lat_s=2
    int size(lat_x * lat_t * lat_spin);
    int MAX_ITER(1e2);
    double TOL(1e-5);
    lattice_vec b(size);
    b.clean1();
    b[0] = 10.0;
    lattice_vec U(size);
    U.clean1();
    double mass(1);
    bool dag(true);
    lattice_dslash Dslash(U, lat_x, lat_t, lat_spin, mass, dag);
    lattice_bicg Bicg(MAX_ITER, TOL, b, Dslash);

    clock_t start = clock();
    Bicg.solve();
    clock_t end = clock();
    std::cout
        << "################"
        << "time cost:"
        << (double)(end - start) / CLOCKS_PER_SEC
        << "s"
        << std::endl;
    return 0;
}
