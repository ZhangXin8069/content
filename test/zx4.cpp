#ifndef _LATTICE_DSLASH_HPP
#define _LATTICE_DSLASH_HPP
#include <complex.h>
#include "lattice_vec.hpp"
#include <mpi.h>
class lattice_dslash
{
public:
    lattice_vec U;
    int lat_x;
    int lat_t;
    int lat_spin;
    double mass;
    bool dag;
    int rank_;
    int size_;

    lattice_dslash(lattice_vec &U, int &lat_x, int &lat_t, int &lat_spin, double &mass, bool &dag) : U(U), lat_x(lat_x), lat_t(lat_t), lat_spin(lat_spin), mass(mass), dag(dag)
    {
    }

    ~lattice_dslash()
    {
        // delete[] send_buf;
        // delete[] recv_buf;
    }
    void dslash(lattice_vec &src, lattice_vec &dest)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        dest.clean();
        double a = 2.0;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> tmp;
        double Half = 0.5;
        double flag = (dag == true) ? -1 : 1;
        int lat_t_ = lat_t / size_;
        for (int x = 0; x < lat_x; x++)
        {
            MPI_Barrier(MPI_COMM_WORLD);

            for (int t = 0; t < lat_t_; t++)
            {
                int t0 = lat_t_ * rank_ + t;
                // mass term
                for (int s = 0; s < lat_spin; s++)
                {
                    dest[(x * lat_t_ + t) * 2 + s] += -(a + mass) * src[(x * lat_t_ + t) * 2 + s];
                }

                // backward x
                int b_x = (x + lat_x - 1) % lat_x;
                tmp = (src[(x * lat_t_ + t) * 2 + 0] + flag * src[(x * lat_t_ + t) * 2 + 1]) * Half * U[(b_x * lat_t + t0) * 2 + 0];
                dest[(b_x * lat_t_ + t) * 2 + 0] += tmp;
                dest[(b_x * lat_t_ + t) * 2 + 1] += flag * tmp;

                // forward x
                int f_x = (x + 1) % lat_x;
                tmp = (src[(x * lat_t_ + t) * 2 + 0] - flag * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t0) * 2 + 0]);
                dest[(f_x * lat_t_ + t) * 2 + 0] += tmp;
                dest[(f_x * lat_t_ + t) * 2 + 1] -= flag * tmp;

                // backward t
                int b_t = (t + lat_t_ - 1) % lat_t_;
                int b_t0 = (t0 + lat_t - 1) % lat_t;
                tmp = (src[(x * lat_t_ + t) * 2 + 0] + flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * U[(x * lat_t + b_t0) * 2 + 1];
                if (t == 0)
                {
                    int b_rank = (rank_ + size_ - 1) % size_;
                    MPI_Request b_send;
                    MPI_Request b_recv;
                    // MPI_Barrier(MPI_COMM_WORLD);
                    MPI_Isend(&tmp, 1, MPI_DOUBLE_COMPLEX, b_rank, 0, MPI_COMM_WORLD, &b_send);
                    MPI_Irecv(&tmp, 1, MPI_DOUBLE_COMPLEX, b_rank, 0, MPI_COMM_WORLD, &b_recv);
                    // MPI_Wait(&b_recv, MPI_STATUS_IGNORE);
                    std::cout << "b_recv:Rank# " << rank_ << "<-Rank# " << b_rank << std::endl;
                    dest[(x * lat_t_ + t) * 2 + 0] += tmp;
                    dest[(x * lat_t_ + t) * 2 + 1] += flag * i * tmp;
                    MPI_Wait(&b_send, MPI_STATUS_IGNORE);
                    std::cout << "b_send:Rank# " << rank_ << "->Rank# " << b_rank << std::endl;
                }
                else
                {
                    dest[(x * lat_t_ + b_t) * 2 + 0] += tmp;
                    dest[(x * lat_t_ + b_t) * 2 + 1] -= flag * i * tmp;
                }

                // forward t
                int f_t = (t + 1) % lat_t_;
                tmp = (src[(x * lat_t_ + t) * 2 + 0] - flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t0) * 2 + 1]);
                if (t == lat_t_ - 1)
                {
                    int f_rank = (rank_ + 1) % size_;
                    MPI_Request f_send;
                    MPI_Request f_recv;
                    MPI_Isend(&tmp, 1, MPI_DOUBLE_COMPLEX, f_rank, 0, MPI_COMM_WORLD, &f_send);
                    MPI_Irecv(&tmp, 1, MPI_DOUBLE_COMPLEX, f_rank, 0, MPI_COMM_WORLD, &f_recv);
                    MPI_Wait(&f_recv, MPI_STATUS_IGNORE);
                    std::cout << "f_recv:Rank# " << rank_ << "<-Rank# " << f_rank << std::endl;
                    dest[(x * lat_t_ + t) * 2 + 0] += tmp;
                    dest[(x * lat_t_ + t) * 2 + 1] -= flag * i * tmp;
                    MPI_Wait(&f_send, MPI_STATUS_IGNORE);
                    std::cout << "f_send:Rank# " << rank_ << "->Rank# " << f_rank << std::endl;
                }
                else
                {
                    dest[(x * lat_t_ + f_t) * 2 + 0] += tmp;
                    dest[(x * lat_t_ + f_t) * 2 + 1] += flag * i * tmp;
                }
            }
        }
    }

    void __dslash(std::complex<double> *src, std::complex<double> *dest)
    {
        int rank_;
        int size_;
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        double a = 2.0;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> tmp;
        double Half = 0.5;
        double flag = (dag == true) ? -1 : 1;
        for (int x = 0; x < lat_x; x++)
        {
            for (int t = lat_t / size_ * rank_; t < lat_t / size_ * (rank_ + 1); t++)
            {

                // mass term
                for (int s = 0; s < lat_spin; s++)
                {
                    dest[(x * lat_t + t) * 2 + s] += -(a + mass) * src[(x * lat_t + t) * 2 + s];
                }

                // backward x
                int b_x = (x + lat_x - 1) % lat_x;
                tmp = (src[(x * lat_t + t) * 2 + 0] + flag * src[(x * lat_t + t) * 2 + 1]) * Half * U[(b_x * lat_t + t) * 2 + 0];
                dest[(b_x * lat_t + t) * 2 + 0] += tmp;
                dest[(b_x * lat_t + t) * 2 + 1] += flag * tmp;

                // forward x
                int f_x = (x + 1) % lat_x;
                tmp = (src[(x * lat_t + t) * 2 + 0] - flag * src[(x * lat_t + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t) * 2 + 0]);
                dest[(f_x * lat_t + t) * 2 + 0] += tmp;
                dest[(f_x * lat_t + t) * 2 + 1] -= flag * tmp;

                // backward t
                int b_t = (t + lat_t - 1) % lat_t;
                tmp = (src[(x * lat_t + t) * 2 + 0] + flag * i * src[(x * lat_t + t) * 2 + 1]) * Half * U[(x * lat_t + b_t) * 2 + 1];
                dest[(x * lat_t + b_t) * 2 + 0] += tmp;
                dest[(x * lat_t + b_t) * 2 + 1] -= flag * i * tmp;

                // forward t
                int f_t = (t + 1) % lat_t;
                tmp = (src[(x * lat_t + t) * 2 + 0] - flag * i * src[(x * lat_t + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t) * 2 + 1]);
                dest[(x * lat_t + f_t) * 2 + 0] += tmp;
                dest[(x * lat_t + f_t) * 2 + 1] += flag * i * tmp;
            }
        }
    }
    void _dslash(lattice_vec &src, lattice_vec &dest)
    {
        dest = src * 2.0;
        dest[0] = 0.5;
    }
};

#endif