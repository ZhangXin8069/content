#include <complex>
#include <iostream>
#include <Eigen/Dense>
#include <ctime>

using namespace std;
using namespace Eigen;

const int MAX_ITER = 1e6; // maximum number of iterations
const double TOL = 1e-12; // tolerance for convergence

template <typename T>
class BiCGStabSolver
{
public:
    BiCGStabSolver(const int &dim4x, const int &dim4t, const int &dim4s, const Vector<T, Dynamic> &U, const Vector<T, Dynamic> &b, const double &mass, const bool &dag) : dim4x(dim4x), dim4t(dim4t), dim4s(dim4s), U(U), b(b), mass(mass), dag(dag) {}

    // solve the linear system Ax = b
    Vector<T, Dynamic> solve()
    {
        int n = dim4x * dim4t * dim4s;
        Vector<T, Dynamic> x = Vector<T, Dynamic>::Zero(n); // initial guess
        // Vector<T, Dynamic> r = b - A * x;
        Vector<T, Dynamic> r = b;
        Vector<T, Dynamic> r_tilde = r;
        Vector<T, Dynamic> p = Vector<T, Dynamic>::Zero(n);
        Vector<T, Dynamic> v = Vector<T, Dynamic>::Zero(n);
        Vector<T, Dynamic> s = Vector<T, Dynamic>::Zero(n);
        Vector<T, Dynamic> t = Vector<T, Dynamic>::Zero(n);

        T rho_1 = T(1);
        T alpha = T(1);
        T omega = T(1);
        T beta = T(1);

        for (int i = 0; i < MAX_ITER; i++)
        {
            // cout << "#" << i << "-Residual: " << r.norm() << endl;
            T rho_2 = r_tilde.dot(r);
            if (rho_2 == T(0))
            {
                cerr << "BiCGStab failed to converge" << endl;
                break;
            }
            if (i > 0)
            {
                beta = (rho_2 / rho_1) * (alpha / omega);
                p = r + beta * (p - omega * v);
            }
            else
            {
                p = r;
            }
            // v = A * p;
            dslash(p, v);
            alpha = rho_2 / r_tilde.dot(v);
            s = r - alpha * v;
            x += alpha * p;
            // t = A * s;
            dslash(s, t);
            omega = t.dot(s) / t.dot(t);
            x += alpha * p + omega * s;
            r = s - omega * t;

            rho_1 = rho_2;
            if (r.norm() < TOL)
            {
                break;
            }
        }
        cout << "#End-Residual: " << r.norm() << endl;
        // dslash(x, p);
        // r = b - p;
        // cout << "#End-Residual: " << r.norm() << endl;
        return x;
    }

private:
    Vector<T, Dynamic> U, b;
    double mass;
    bool dag;
    int dim4x;
    int dim4t;
    int dim4s;
    void dslash(const Vector<T, Dynamic> &src, Vector<T, Dynamic> &dest)
    {

        dest = VectorXcd::Zero(dest.size());
        const double a = 2.0;
        const complex<double> i(0.0, 1.0);
        complex<double> tmp;
        const double Half = 0.5;
        double flag = (dag == true) ? -1 : 1;
        for (int x = 0; x < dim4x; x++)
            for (int t = 0; t < dim4t; t++)
            {

                // mass term
                for (int s = 0; s < dim4s; s++)
                {
                    dest((x * dim4t + t) * 2 + s) += -(a + mass) * src((x * dim4t + t) * 2 + s);
                }

                // backward x
                int b_x = (x + dim4x - 1) % dim4x;
                tmp = (src((x * dim4t + t) * 2 + 0) + flag * src((x * dim4t + t) * 2 + 1)) * Half * U((b_x * dim4t + t) * 2 + 0);
                dest((b_x * dim4t + t) * 2 + 0) += tmp;
                dest((b_x * dim4t + t) * 2 + 1) += flag * tmp;

                // forward x
                int f_x = (x + 1) % dim4x;
                tmp = (src((x * dim4t + t) * 2 + 0) - flag * src((x * dim4t + t) * 2 + 1)) * Half * conj(U((x * dim4t + t) * 2 + 0));
                dest((f_x * dim4t + t) * 2 + 0) += tmp;
                dest((f_x * dim4t + t) * 2 + 1) -= flag * tmp;

                // backward t
                int b_t = (t + dim4t - 1) % dim4t;
                tmp = (src((x * dim4t + t) * 2 + 0) + flag * i * src((x * dim4t + t) * 2 + 1)) * Half * U((x * dim4t + b_t) * 2 + 1);
                dest((x * dim4t + b_t) * 2 + 0) += tmp;
                dest((x * dim4t + b_t) * 2 + 1) -= flag * i * tmp;

                // forward t
                int f_t = (t + 1) % dim4t;
                tmp = (src((x * dim4t + t) * 2 + 0) - flag * i * src((x * dim4t + t) * 2 + 1)) * Half * conj(U((x * dim4t + t) * 2 + 1));
                dest((x * dim4t + f_t) * 2 + 0) += tmp;
                dest((x * dim4t + f_t) * 2 + 1) += flag * i * tmp;
            }
    }
};

int main()
{

    double mass(1);
    bool dag(true);
    int dim4x(40);
    int dim4t(40);
    int dim4s(2);
    const int nu4dim(3200);
    Vector<complex<double>, nu4dim> U;
    Vector<complex<double>, nu4dim> b;
    U = VectorXcd::Ones(nu4dim);
    // U = VectorXcd::Random(nu4dim);
    b = VectorXcd::Ones(nu4dim);
    b(nu4dim-1)=10;
    // b = VectorXcd::Random(nu4dim);
    clock_t start = clock();
    BiCGStabSolver<complex<double>> solver(dim4x, dim4t, dim4s, U, b, mass, dag);
    Vector<complex<double>, nu4dim> x = solver.solve();
    cout << "Solution:" << x << endl;
    clock_t end = clock();
    cout
        // << "rank:"
        // << rank
        << "################"
        << "time cost:"
        << (double)(end - start) / CLOCKS_PER_SEC
        << "s"
        << endl;
    return 0;
}