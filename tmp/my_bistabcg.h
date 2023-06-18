#ifndef _MY_BISTABCG_H
#define _MY_BISTABCG_H
#include "my_vector.h"
#include "my_dslash.h"
#include <iostream>
class BistabCG
{
public:
    BistabCG(ComplexVector &U,
             ComplexVector &b,
             int &lattice_x,
             int &lattice_t,
             int &lattice_spin,
             double &mass,
             bool &dag,
             int &MAX_ITER,
             const double &TOL) : U_(U),
                                  b_(b),
                                  lattice_x_(lattice_x),
                                  lattice_t_(lattice_t),
                                  lattice_spin_(lattice_spin),
                                  mass_(mass),
                                  dag_(dag),
                                  MAX_ITER_(MAX_ITER),
                                  TOL_(TOL)
    {
        MPI_Comm_size(MPI_COMM_WORLD, &node_size_);
        MPI_Comm_rank(MPI_COMM_WORLD, &node_rank_);
        lattice_x_ = lattice_x / node_size_;
    }

    ~BistabCG()
    {
    }

    // solve the linear system Ax = b_
    void solve()
    {
        int size(lattice_x_ * lattice_t_ * lattice_spin_);
        std::complex<double> rho_prev(1.0, 0.0), rho(0.0, 0.0), alpha(1.0, 0.0), omega(1.0, 0.0), beta(0.0, 0.0);
        ComplexVector x(size), r(size), r_tilde(size), p(size), v(size), s(size), t(size);

        for (int x = 0; x < lattice_x_; x++)
        {
            int x0 = lattice_x_ * node_rank_ + x;
            for (int t = 0; t < lattice_t_; t++)
            {
                for (int s = 0; s < lattice_spin_; s++)
                {
                    r[(x * lattice_t_ + t) * 2 + s] = b_[(x0 * lattice_t_ + t) * 2 + s];
                }
            }
        }

        r_tilde = r;
        // x.rand(); // initial guess
        // Dslash_.dslash(x, r_tilde);
        // // ComplexVector r = b_ - A * x;
        // r = r - r_tilde;
        // if x=0;r_tilde = r0 = b_;
        Dslash Dslash_(U_, lattice_x_, lattice_t_, lattice_spin_, mass_, dag_);

        for (int i = 0; i < MAX_ITER_; i++)
        {
            rho = r_tilde.dotX(r);

            beta = (rho / rho_prev) * (alpha / omega);
            p = r + (p - v * omega) * beta;

            // v = A * p;
            Dslash_.dslash(p, v);

            alpha = rho / r_tilde.dotX(v);

            s = r - v * alpha;

            // t = A * s;
            Dslash_.dslash(s, t);

            omega = t.dotX(s) / t.dotX(t);

            x = x + p * alpha + s * omega;

            r = s - t * omega;

            if (r.norm2X().real() < TOL_ || i == MAX_ITER_ - 1)
            {
                std::cout << "##loop "
                          << i
                          << "##Residual:"
                          << r.norm2X().real()
                          << std::endl;
                if (node_rank_ == 0)
                {
                    std::cout << x;
                }
                break;
            }
            rho_prev = rho;
        }
    }

private:
    ComplexVector U_, b_;
    int MAX_ITER_, lattice_x_, lattice_t_, lattice_spin_, node_size_, node_rank_;
    double TOL_, mass_;
    bool dag_;
};

#endif
