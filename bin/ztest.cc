#include "../include/zudaX.h"
int main(int argc, char **argv)
{
    double start, end;
    int lat_x(16), lat_y(16), lat_z(16), lat_t(32), lat_s(4), lat_c(3);
    int num_x(1), num_y(1), num_z(2), num_t(4);
    int MAX_ITER(1e6);
    double TOL(1e-12);
    LatticeGauge U(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    LatticeFermi b(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    LatticeFermi x(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    U.assign_random(000);
    U.assign_unit();
    b.assign_random(111);
    b.assign_zero();
    b(0, 0, 0, 0, 0, 0) = 1.0;
    x.assign_random(222);
    x.assign_zero();
    MPI_Init(&argc, &argv);
    LatticeGauge block_U = U.block(num_x, num_y, num_z, num_t);
    LatticeFermi block_b = b.block(num_x, num_y, num_z, num_t);
    LatticeFermi block_x = x.block(num_x, num_y, num_z, num_t);
    start = MPI_Wtime();
    cg(block_U, block_b, block_x, num_x, num_y, num_z, num_t, MAX_ITER, TOL, true);
    // dslash(block_U, block_b, block_x, false);
    end = MPI_Wtime();
    LatticeFermi reback_x = block_x.reback(num_x, num_y, num_z, num_t);
    x = reback_x - b;
    x.print(0, 0, 0, 0, 0, 0);
    x.print(0, 0, 0, 0, 0, 1);
    // int rank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // if(rank==0){
    // x.print();
    // }
    std::cout
        << "######time cost:" << end - start << "s ######" << std::endl;
    MPI_Finalize();
    return 0;
}