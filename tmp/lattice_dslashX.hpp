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
    int lat_t_;
    int lat_spin;
    double mass;
    bool dag;
    int rank_;
    int size_;
    lattice_dslash(lattice_vec &U, int &lat_x, int &lat_t, int &lat_spin, double &mass, bool &dag) : U(U), lat_x(lat_x), lat_t(lat_t), lat_spin(lat_spin), mass(mass), dag(dag)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
        i = (0.0, 1.0);
        a = 2.0;
        Half = 0.5;
        flag = (dag == true) ? -1 : 1;
        lat_t_ = lat_t / size_;
        b_rank = (rank_ + size_ - 1) % size_;
        f_rank = (rank_ + 1) % size_;
    }

    ~lattice_dslash()
    {
    }
    void _dslash(lattice_vec &src, lattice_vec &dest)
    {
        if (size_ == 1)
        {
            dslash_old(src, dest);
        }
        else
        {
            // dslash_test(src, dest);
            dest.clean();
            dslash_x(src, dest, 0);
            dslash_b_begin(src, 0);
            dslash_x(src, dest, lat_t_ - 1);
            dslash_f_begin(src, lat_t_ - 1);
            for (int t = 1; t < lat_t_ - 1; t++)
            {
                dslash_prev(src, dest, t);
            }
            std::cout << "norm done:Rank# " << rank_ << "######" << std::endl;
            dslash_b_end(dest, 0);
            dslash_f_end(dest, lat_t_ - 1);
        }
    }
    void dslash(lattice_vec &src, lattice_vec &dest)
    {
        dest = src * 2.0;
        dest[0] = 0.5;
    }

private:
    double a;
    std::complex<double> i;
    std::complex<double> tmp;
    double Half;
    double flag;
    int b_rank;
    int f_rank;
    void dslash_prev(lattice_vec &src, lattice_vec &dest, int t)
    {
        int t0 = lat_t_ * rank_ + t;
        for (int x = 0; x < lat_x; x++)
        {
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
            dest[(x * lat_t_ + b_t) * 2 + 0] += tmp;
            dest[(x * lat_t_ + b_t) * 2 + 1] -= flag * i * tmp;

            // forward t
            int f_t = (t + 1) % lat_t_;
            tmp = (src[(x * lat_t_ + t) * 2 + 0] - flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t0) * 2 + 1]);
            dest[(x * lat_t_ + f_t) * 2 + 0] += tmp;
            dest[(x * lat_t_ + f_t) * 2 + 1] += flag * i * tmp;
        }
    }

    void dslash_x(lattice_vec &src, lattice_vec &dest, int t)
    {
        int t0 = lat_t_ * rank_ + t;
        for (int x = 0; x < lat_x; x++)
        {
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
        }
    }
    void dslash_b_begin(lattice_vec &src, int t)
    {
        std::complex<double> b_send_vec[lat_x];
        MPI_Request b_send_req;
        int t0 = lat_t_ * rank_ + t;
        int b_t = (t + lat_t_ - 1) % lat_t_;
        int b_t0 = (t0 + lat_t - 1) % lat_t;
        for (int x = 0; x < lat_x; x++)
        {
            b_send_vec[x] = (src[(x * lat_t_ + t) * 2 + 0] + flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * U[(x * lat_t + b_t0) * 2 + 1];
        }
        MPI_Isend(&b_send_vec, lat_x, MPI_DOUBLE_COMPLEX, b_rank, b_rank, MPI_COMM_WORLD, &b_send_req);
        std::cout << "b_send:Rank# " << rank_ << "->Rank# " << b_rank << std::endl;
    }
    void dslash_f_begin(lattice_vec &src, int t)
    {
        std::complex<double> f_send_vec[lat_x];
        MPI_Request f_send_req;
        int t0 = lat_t_ * rank_ + t;
        int f_t = (t + 1) % lat_t_;
        for (int x = 0; x < lat_x; x++)
        {
            f_send_vec[x] = (src[(x * lat_t_ + t) * 2 + 0] - flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t0) * 2 + 1]);
        }
        MPI_Isend(&f_send_vec, lat_x, MPI_DOUBLE_COMPLEX, f_rank, f_rank, MPI_COMM_WORLD, &f_send_req);
        std::cout << "f_send:Rank# " << rank_ << "->Rank# " << f_rank << std::endl;
    }
    void dslash_b_end(lattice_vec &dest, int t)
    {
        std::complex<double> b_recv_vec[lat_x];
        MPI_Request b_recv_req;
        MPI_Status b_recv_stat;
        MPI_Irecv(&b_recv_vec, lat_x, MPI_DOUBLE_COMPLEX, b_rank, rank_, MPI_COMM_WORLD, &b_recv_req);
        MPI_Wait(&b_recv_req, &b_recv_stat);
        std::cout << "b_recv:Rank# " << rank_ << "<-Rank# " << b_rank << std::endl;
        for (int x = 0; x < lat_x; x++)
        {
            dest[(x * lat_t_ + t) * 2 + 0] += b_recv_vec[x];
            dest[(x * lat_t_ + t) * 2 + 1] += flag * i * b_recv_vec[x];
        }
    }
    void dslash_f_end(lattice_vec &dest, int t)
    {
        std::complex<double> f_recv_vec[lat_x];
        MPI_Request f_recv_req;
        MPI_Status f_recv_stat;
        MPI_Irecv(&f_recv_vec, lat_x, MPI_DOUBLE_COMPLEX, f_rank, rank_, MPI_COMM_WORLD, &f_recv_req);
        MPI_Wait(&f_recv_req, &f_recv_stat);
        std::cout << "f_recv:Rank# " << rank_ << "<-Rank# " << f_rank << std::endl;
        for (int x = 0; x < lat_x; x++)
        {
            dest[(x * lat_t_ + t) * 2 + 0] += f_recv_vec[x];
            dest[(x * lat_t_ + t) * 2 + 1] -= flag * i * f_recv_vec[x];
        }
    }
    void dslash_test(lattice_vec &src, lattice_vec &dest)
    {
        dest = src * 2.0;
        dest[0] = 0.5;
    }
    void dslash_old(lattice_vec &src, lattice_vec &dest)
    {
        dest.clean();
        double a = 2.0;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> tmp;
        double Half = 0.5;
        double flag = (dag == true) ? -1 : 1;
        for (int x = 0; x < lat_x; x++)
            for (int t = 0; t < lat_t; t++)
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
};

#endif