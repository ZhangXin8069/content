#include "lattice_vec.hpp"
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
        i = (0.0, 1.0);
        a = 2.0;
        Half = 0.5;
        flag = (dag == true) ? -1 : 1;
        lat_t_ = lat_t / size_;
        b_rank = (rank_ + size_ - 1) % size_;
        f_rank = (rank_ + 1) % size_;
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
        lattice_vec r = *b_;
        lattice_vec r_tilde(size);
        lattice_vec p(size);
        lattice_vec v(size);
        lattice_vec s(size);
        lattice_vec t(size);
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
            rho = r_tilde.dot(r);
            beta = (rho / rho_prev) * (alpha / omega);
            p = r + (p - v * omega) * beta;
            // v = A * p;
            dslash(p, v);
            alpha = rho / r_tilde.dot(v);
            s = r - v * alpha;
            // t = A * s;
            dslash(s, t);
            omega = t.dot(s) / t.dot(t);
            x = x + p * alpha + s * omega;
            r = s - t * omega;
            if (r.norm() < TOL || i == MAX_ITER - 1)
            {
                std::cout << "##loop "
                          << i
                          << "##Residual:"
                          << r.norm()
                          << std::endl;
                break;
            }
            rho_prev = rho;
        }
        x.print();
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

    void dslash(lattice_vec &src, lattice_vec &dest)
    {
        dest = src * 2.0 + 0.5;
    }
};

int main()
{
    int lat_x(16);
    int lat_t(16);
    int lat_spin(2); // const lat_spin=2
    int size(lat_x * lat_t * lat_spin);
    int MAX_ITER(1e5);
    double TOL(1e-5);
    lattice_vec b(size);
    lattice_vec U(size);
    b.clean1();
    b[0] = 10;
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
    return 0;
}