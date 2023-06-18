#ifndef _LATTICE_BICG_HPP
#define _LATTICE_BICG_HPP
#include "lattice_vec.hpp"
#include "lattice_dslash.hpp"
class lattice_bicg
{
public:
    lattice_bicg(const int &MAX_ITER, const double &TOL, lattice_vec &b, lattice_dslash &Dslash) : MAX_ITER(MAX_ITER), TOL(TOL), Dslash(Dslash)
    {
        this->size = b.size;
        this->b = &b;
    }

    ~lattice_bicg()
    {
    }

    // solve the linear system Ax = b
    void solve()
    {
        std::complex<double> rho_prev(1.0, 0.0);
        std::complex<double> rho(0.0, 0.0);
        std::complex<double> alpha(1.0, 0.0);
        std::complex<double> omega(1.0, 0.0);
        std::complex<double> beta(0.0, 0.0);
        lattice_vec x(size);
        lattice_vec r = *b;
        lattice_vec r_tilde(size);
        lattice_vec p(size);
        lattice_vec v(size);
        lattice_vec s(size);
        lattice_vec t(size);
        // x.rand(); // initial guess
        // Dslash.dslash(x, r_tilde);
        // // lattice_vec r = b - A * x;
        // r = r - r_tilde;
        x.clean();//if x=0;r_tilde=r0=b;
        r_tilde = r;
        p.clean();
        v.clean();
        s.clean();
        t.clean();
        for (int i = 0; i < MAX_ITER; i++)
        {
            std::cout << "#" << i << "-Residual: " << r.norm() << std::endl;
            rho = r_tilde.dot(r);
            beta = (rho / rho_prev) * (alpha / omega);
            p = r + (p - v * omega) * beta;
            // v = A * p;
            Dslash.dslash(p, v);
            alpha = rho / r_tilde.dot(v);
            s = r - v * alpha;
            // t = A * s;
            Dslash.dslash(s, t);
            omega = t.dot(s) / t.dot(t);
            x = x + p * alpha + s * omega;
            r = s - t * omega;
            if (r.norm() < TOL)
            {
                break;
            }
            rho_prev = rho;
        }
        std::cout << "#End-Residual: " << r.norm() << std::endl;
        x.print();
        // Dslash.dslash(x, p);
        // r = b - p;
        // std::cout << "#End-Residual: " << r.norm() << std::endl;
    }
private:
    int MAX_ITER;
    double TOL;
    lattice_vec *b;
    int size;
    lattice_dslash Dslash;
};
#endif