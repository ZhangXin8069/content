#include "lattice_vec.hpp"
#include <mpi.h>
#include <ctime>
#include <iostream>
#include <complex.h>
class lattice_bicg
{
public:
    lattice_bicg(lattice_vec &U,
                 lattice_vec &b,
                 int &lat_x,
                 int &lat_t,
                 int &lat_spin,
                 double &mass,
                 bool &dag,
                 int &MAX_ITER,
                 const double &TOL) : U_(&U),
                                      b_(&b), lat_x(lat_x),
                                      lat_t(lat_t),
                                      lat_spin(lat_spin),
                                      mass(mass),
                                      dag(dag),
                                      MAX_ITER(MAX_ITER),
                                      TOL(TOL),
                                      size(lat_x * lat_t * lat_spin)
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
        b_send_vec = new std::complex<double>[lat_x];
        b_recv_vec = new std::complex<double>[lat_x];
        f_send_vec = new std::complex<double>[lat_x];
        f_recv_vec = new std::complex<double>[lat_x];
    }

    ~lattice_bicg()
    {
        delete[] b_send_vec;
        delete[] b_recv_vec;
        delete[] f_send_vec;
        delete[] f_recv_vec;
    }

    // solve the linear system Ax = b
    void solve()
    {
        size = size / size_;
        std::complex<double> rho_prev(1.0, 0.0);
        std::complex<double> rho(0.0, 0.0);
        std::complex<double> tmp(0.0, 0.0);
        std::complex<double> tmpX(0.0, 0.0);
        double tmpXX;
        std::complex<double> alpha(1.0, 0.0);
        std::complex<double> omega(1.0, 0.0);
        std::complex<double> beta(0.0, 0.0);
        lattice_vec x(size);
        lattice_vec r(size);
        lattice_vec r_tilde(size);
        lattice_vec p(size);
        lattice_vec v(size);
        lattice_vec s(size);
        lattice_vec t(size);
        for (int x = 0; x < lat_x; x++)
        {
            for (int t = 0; t < lat_t / size_; t++)
            {
                int t0 = lat_t / size_ * rank_ + t;
                for (int s = 0; s < lat_spin; s++)
                {
                    r[(x * lat_t / size_ + t) * 2 + s] = (*b_)[(x * lat_t + t0) * 2 + s];
                }
            }
        }
        // x.rand(); // initial guess
        // dslash(x, r_tilde);
        // // lattice_vec r = b - A * x;
        // r = r - r_tilde;
        x.clean(); // if x=0;r_tilde=r0=b;
        r_tilde = r;
        p.clean();
        v.clean();
        s.clean();
        t.clean();
        for (int i = 0; i < MAX_ITER; i++)
        {
            dotX(rho, r_tilde, r);

            beta = (rho / rho_prev) * (alpha / omega);
            p = r + (p - v * omega) * beta;
            if (rank_ == 0)
            {
                // std::cout << "##loop "
                //           << i
                //           << "##Test:"
                //           << tmp
                //           << std::endl;
                p.print();
            }
            // v = A * p;
            dslash(p, v);

            dotX(tmp, r_tilde, v);
            alpha = rho / tmp;

            s = r - v * alpha;

            // t = A * s;
            dslash(s, t);

            dotX(tmp, t, s);
            dotX(tmpX, t, t);
            omega = tmp / tmpX;

            x = x + p * alpha + s * omega;
            r = s - t * omega;

            normX(tmpXX, r);
            if (tmpXX < TOL || i == MAX_ITER - 1)
            {
                std::cout << "##loop "
                          << i
                          << "##Residual:"
                          << tmpXX
                          << std::endl;
                break;
            }

            rho_prev = rho;
        }
        if (rank_ == 0)
        {
            x.print();
        }
    }

private:
    lattice_vec *U_;
    lattice_vec *b_;
    int MAX_ITER;
    double TOL;
    int lat_x;
    int lat_t;
    int lat_t_;
    int lat_spin;
    double mass;
    bool dag;
    int size;
    int rank_;
    int size_;
    double a;
    std::complex<double> i;
    std::complex<double> tmp;
    double Half;
    double flag;
    int b_rank;
    int f_rank;
    std::complex<double> *b_send_vec;
    std::complex<double> *b_recv_vec;
    std::complex<double> *f_send_vec;
    std::complex<double> *f_recv_vec;
    void dslash(lattice_vec &src, lattice_vec &dest)
    {
        if (size_ == 1)
        {
            dslash_old(src, dest);
        }
        else if (1)
        {
            dslash_test(src, dest);
        }
        else
        {
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
    void dslash_test(lattice_vec &src, lattice_vec &dest)
    {
        dest = src * 2.0 + 0.5;
    }
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
            tmp = (src[(x * lat_t_ + t) * 2 + 0] + flag * src[(x * lat_t_ + t) * 2 + 1]) * Half * (*U_)[(b_x * lat_t + t0) * 2 + 0];
            dest[(b_x * lat_t_ + t) * 2 + 0] += tmp;
            dest[(b_x * lat_t_ + t) * 2 + 1] += flag * tmp;

            // forward x
            int f_x = (x + 1) % lat_x;
            tmp = (src[(x * lat_t_ + t) * 2 + 0] - flag * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj((*U_)[(x * lat_t + t0) * 2 + 0]);
            dest[(f_x * lat_t_ + t) * 2 + 0] += tmp;
            dest[(f_x * lat_t_ + t) * 2 + 1] -= flag * tmp;

            // backward t
            int b_t = (t + lat_t_ - 1) % lat_t_;
            int b_t0 = (t0 + lat_t - 1) % lat_t;
            tmp = (src[(x * lat_t_ + t) * 2 + 0] + flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * (*U_)[(x * lat_t + b_t0) * 2 + 1];
            dest[(x * lat_t_ + b_t) * 2 + 0] += tmp;
            dest[(x * lat_t_ + b_t) * 2 + 1] -= flag * i * tmp;

            // forward t
            int f_t = (t + 1) % lat_t_;
            tmp = (src[(x * lat_t_ + t) * 2 + 0] - flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj((*U_)[(x * lat_t + t0) * 2 + 1]);
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
            tmp = (src[(x * lat_t_ + t) * 2 + 0] + flag * src[(x * lat_t_ + t) * 2 + 1]) * Half * (*U_)[(b_x * lat_t + t0) * 2 + 0];
            dest[(b_x * lat_t_ + t) * 2 + 0] += tmp;
            dest[(b_x * lat_t_ + t) * 2 + 1] += flag * tmp;

            // forward x
            int f_x = (x + 1) % lat_x;
            tmp = (src[(x * lat_t_ + t) * 2 + 0] - flag * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj((*U_)[(x * lat_t + t0) * 2 + 0]);
            dest[(f_x * lat_t_ + t) * 2 + 0] += tmp;
            dest[(f_x * lat_t_ + t) * 2 + 1] -= flag * tmp;
        }
    }
    void dslash_b_begin(lattice_vec &src, int t)
    {
        // std::complex<double> b_send_vec[lat_x];
        //????
        MPI_Request b_send_req;
        int t0 = lat_t_ * rank_ + t;
        int b_t = (t + lat_t_ - 1) % lat_t_;
        int b_t0 = (t0 + lat_t - 1) % lat_t;
        for (int x = 0; x < lat_x; x++)
        {
            b_send_vec[x] = (src[(x * lat_t_ + t) * 2 + 0] + flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * (*U_)[(x * lat_t + b_t0) * 2 + 1];
        }
        MPI_Isend(&b_send_vec, lat_x, MPI_DOUBLE_COMPLEX, b_rank, b_rank, MPI_COMM_WORLD, &b_send_req);
        std::cout << "b_send:Rank# " << rank_ << "->Rank# " << b_rank << std::endl;
    }
    void dslash_f_begin(lattice_vec &src, int t)
    {
        // std::complex<double> f_send_vec[lat_x];
        MPI_Request f_send_req;
        int t0 = lat_t_ * rank_ + t;
        int f_t = (t + 1) % lat_t_;
        for (int x = 0; x < lat_x; x++)
        {
            f_send_vec[x] = (src[(x * lat_t_ + t) * 2 + 0] - flag * i * src[(x * lat_t_ + t) * 2 + 1]) * Half * conj((*U_)[(x * lat_t + t0) * 2 + 1]);
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
                tmp = (src[(x * lat_t + t) * 2 + 0] + flag * src[(x * lat_t + t) * 2 + 1]) * Half * (*U_)[(b_x * lat_t + t) * 2 + 0];
                dest[(b_x * lat_t + t) * 2 + 0] += tmp;
                dest[(b_x * lat_t + t) * 2 + 1] += flag * tmp;

                // forward x
                int f_x = (x + 1) % lat_x;
                tmp = (src[(x * lat_t + t) * 2 + 0] - flag * src[(x * lat_t + t) * 2 + 1]) * Half * conj((*U_)[(x * lat_t + t) * 2 + 0]);
                dest[(f_x * lat_t + t) * 2 + 0] += tmp;
                dest[(f_x * lat_t + t) * 2 + 1] -= flag * tmp;

                // backward t
                int b_t = (t + lat_t - 1) % lat_t;
                tmp = (src[(x * lat_t + t) * 2 + 0] + flag * i * src[(x * lat_t + t) * 2 + 1]) * Half * (*U_)[(x * lat_t + b_t) * 2 + 1];
                dest[(x * lat_t + b_t) * 2 + 0] += tmp;
                dest[(x * lat_t + b_t) * 2 + 1] -= flag * i * tmp;

                // forward t
                int f_t = (t + 1) % lat_t;
                tmp = (src[(x * lat_t + t) * 2 + 0] - flag * i * src[(x * lat_t + t) * 2 + 1]) * Half * conj((*U_)[(x * lat_t + t) * 2 + 1]);
                dest[(x * lat_t + f_t) * 2 + 0] += tmp;
                dest[(x * lat_t + f_t) * 2 + 1] += flag * i * tmp;
            }
    }
    void dotX(std::complex<double> &tmp, lattice_vec &t, lattice_vec &s)
    {
        tmp = t.dot(s);
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmp, &tmp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    }
    void normX(double &tmpXX, lattice_vec &r)
    {
        tmpXX = r.norm();
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmpXX, &tmpXX, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
    }
};

int main()
{
    MPI_Init(NULL, NULL);
    int lat_x(100);
    int lat_t(100);
    int lat_spin(2); // const lat_spin=2
    int size(lat_x * lat_t * lat_spin);
    int MAX_ITER(1e1);
    double TOL(1e-5);
    lattice_vec b(size);
    b.clean1();
    b[0] = 10.0;
    lattice_vec U(size);
    U.clean1();
    double mass(1);
    bool dag(true);
    clock_t start = clock();
    lattice_bicg Bicg(U, b, lat_x, lat_t, lat_spin, mass, dag, MAX_ITER, TOL);
    Bicg.solve();
    clock_t end = clock();
    std::cout
        << "################"
        << "time cost:"
        << (double)(end - start) / CLOCKS_PER_SEC
        << "s"
        << std::endl;
    MPI_Finalize();
    return 0;
}