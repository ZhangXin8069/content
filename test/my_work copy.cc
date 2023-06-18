#include "my_bistabcg.h"
int main()
{
    MPI_Init(NULL, NULL);
    int lattice_x(60), lattice_t(60), lattice_spin(2); // const lattice_spin=2
    int size(lattice_x * lattice_t * lattice_spin);
    int MAX_ITER(1e3);
    double TOL(1e-5), start, end;
    ComplexVector b(size);
    ComplexVector U(size);
    // b.clean_1();
    // // b[0] = 10;
    // U.clean_1();
    b.clean_rand();
    U.clean_rand();
    double mass(1);
    bool dag(true);
    start = MPI_Wtime();
    BistabCG bicgstaCG(U, b, lattice_x, lattice_t, lattice_spin, mass, dag, MAX_ITER, TOL);
    bicgstaCG.solve();
    end = MPI_Wtime();
    std::cout
        << "########"
        << "time cost:"
        << end - start
        << "s"
        << std::endl;
    MPI_Finalize();

    return 0;
}