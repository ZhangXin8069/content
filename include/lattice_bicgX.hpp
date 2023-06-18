#ifndef _LATTICE_BICG_HPP
#define _LATTICE_BICG_HPP
#include "lattice_vec.hpp"
#include "lattice_dslashX.hpp"
class lattice_bicg
{
public:
    lattice_bicg(const int &MAX_ITER, const double &TOL, lattice_vec &b, lattice_dslash &Dslash) : MAX_ITER(MAX_ITER), TOL(TOL), Dslash(Dslash)
    {
        this->size = b.size;
        this->b = &b;
        MPI_Comm_size(MPI_COMM_WORLD, &size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &rank_);
    }

    ~lattice_bicg()
    {
    }

    // solve the linear system Ax = b
    void solve()
    {
        lattice_vec b0 = *b;
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
        for (int x = 0; x < Dslash.lat_x; x++)
        {
            for (int t = 0; t < Dslash.lat_t / size_; t++)
            {
                int t0 = Dslash.lat_t / size_ * rank_ + t;
                for (int s = 0; s < Dslash.lat_spin; s++)
                {
                    r[(x * Dslash.lat_t / size_ + t) * 2 + s] = b0[(x * Dslash.lat_t + t0) * 2 + s];
                }
            }
        }
        // x.rand(); // initial guess
        // Dslash.dslash(x, r_tilde);
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

            // v = A * p;
            Dslash.dslash(p, v);

            dotX(tmp, r_tilde, v);
            alpha = rho / tmp;

            s = r - v * alpha;

            // t = A * s;
            Dslash.dslash(s, t);

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
        // x.delete_();
        // // r.delete_();//r==b
        // r_tilde.delete_();
        // p.delete_();
        // v.delete_();
        // s.delete_();
        // t.delete_();
    }

private:
    int MAX_ITER;
    double TOL;
    lattice_vec *b;
    int size;
    lattice_dslash Dslash;
    int rank_;
    int size_;
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
#endif