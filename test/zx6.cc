#include <complex>
#include <vector>
#include <math.h>
#include <ctime>
#include <iostream>
#include "Eigen/Core"
using namespace std;
using namespace Eigen;

class lattice_fermi
{
public:
    int lat_x, lat_t, lat_spin;
    int size;
    vector<complex<double>> lattice_vec;
    lattice_fermi(int lat_x, int lat_t, int lat_spin) : lat_x(lat_x), lat_t(lat_t), lat_spin(lat_spin), lattice_vec(lat_x * lat_t * lat_spin) { size = lattice_vec.size(); }
    lattice_fermi();
    void clean()
    {
        for (int i = 0; i < size; i++)
            lattice_vec[i] = 0;
    }
    complex<double> &operator[](int i) { return lattice_vec[i]; }
    lattice_fermi &operator=(lattice_fermi &a)
    {
        for (int i = 0; i < size; i++)
            this->lattice_vec[i] = a.lattice_vec[i];
        return *this;
    }
    lattice_fermi &operator-(const lattice_fermi &a)
    {
        for (int i = 0; i < size; i++)
            this->lattice_vec[i] = this->lattice_vec[i] - a.lattice_vec[i];
        return *this;
    }
    lattice_fermi &operator+(const lattice_fermi &a)
    {
        for (int i = 0; i < size; i++)
            this->lattice_vec[i] = this->lattice_vec[i] + a.lattice_vec[i];
        return *this;
    }
};
class lattice_gauge
{
public:
    int lat_x, lat_t, lat_d;
    int size;
    vector<complex<double>> lattice_vec_c;
    lattice_gauge(int lat_x, int lat_t, int lat_d) : lat_x(lat_x), lat_t(lat_t), lat_d(lat_d), lattice_vec_c(lat_x * lat_t * lat_d) { size = lattice_vec_c.size(); }
    lattice_gauge();
    complex<double> &operator[](int i) { return lattice_vec_c[i]; }
    lattice_gauge &operator=(lattice_gauge a)
    {
        for (int i = 0; i < size; i++)
            this->lattice_vec_c[i] = a.lattice_vec_c[i];
        return *this;
    }
};
class lattice_propagator
{
public:
    int lat_x, lat_t, lat_spin;
    int size;
    vector<complex<double>> lattice_vec_c;
    lattice_propagator(int lat_x, int lat_t, int lat_spin) : lat_x(lat_x), lat_t(lat_t), lat_spin(lat_spin), lattice_vec_c(lat_x * lat_t * lat_spin * lat_spin) { size = lattice_vec_c.size(); }
    lattice_propagator();
    void clean()
    {
        for (int i = 0; i < size; i++)
            lattice_vec_c[i] = 0;
    }
    complex<double> &operator[](int i) { return lattice_vec_c[i]; }
    lattice_propagator &operator=(lattice_propagator a)
    {
        for (int i = 0; i < size; i++)
            this->lattice_vec_c[i] = a.lattice_vec_c[i];
        return *this;
    }
};
class Gamma
{
public:
    int lat_spin;
    int size;
    vector<complex<double>> lattice_vec_c;
    Gamma(int lat_spin) : lat_spin(lat_spin), lattice_vec_c(lat_spin * lat_spin) { size = lattice_vec_c.size(); }
    Gamma();
    void clean()
    {
        for (int i = 0; i < size; i++)
            lattice_vec_c[i] = 0;
    }
    complex<double> &operator[](int i) { return lattice_vec_c[i]; }
    Gamma &operator=(Gamma a)
    {
        for (int i = 0; i < size; i++)
            this->lattice_vec_c[i] = a.lattice_vec_c[i];
        return *this;
    }
};
lattice_propagator operator*(Gamma G, lattice_propagator prop)
{
    lattice_propagator prop1(prop);
    prop1.clean();
    for (int x = 0; x < prop.lat_x; x++)
        for (int t = 0; t < prop.lat_t; t++)
            for (int s1 = 0; s1 < prop.lat_spin; s1++)
                for (int s2 = 0; s2 < prop.lat_spin; s2++)
                {
                    for (int i = 0; i < prop.lat_spin; i++)
                        prop1[x * prop.lat_t * prop.lat_spin * prop.lat_spin + t * prop.lat_spin * prop.lat_spin + s1 * prop.lat_spin + s2] += G[s1 * G.lat_spin + i] * prop[x * prop.lat_t * prop.lat_spin * prop.lat_spin + t * prop.lat_spin * prop.lat_spin + i * prop.lat_spin + s2];
                }
    return prop1;
}
lattice_propagator operator*(lattice_propagator prop, Gamma G)
{
    lattice_propagator prop1(prop);
    prop1.clean();
    for (int x = 0; x < prop.lat_x; x++)
        for (int t = 0; t < prop.lat_t; t++)
            for (int s1 = 0; s1 < prop.lat_spin; s1++)
                for (int s2 = 0; s2 < prop.lat_spin; s2++)
                {
                    for (int i = 0; i < prop.lat_spin; i++)
                        prop1[x * prop.lat_t * prop.lat_spin * prop.lat_spin + t * prop.lat_spin * prop.lat_spin + s1 * prop.lat_spin + s2] += G[i * G.lat_spin + s2] * prop[x * prop.lat_t * prop.lat_spin * prop.lat_spin + t * prop.lat_spin * prop.lat_spin + s1 * prop.lat_spin + i];
                }
    return prop1;
}
double norm_2(lattice_fermi s)
{
    complex<double> s1(0.0, 0.0);
    for (int i = 0; i < s.size; i++)
    {
        s1 += s[i] * conj(s[i]);
    }
    return s1.real();
};
double norm_2(lattice_propagator f)
{
    complex<double> f1(0.0, 0.0);
    for (int i = 0; i < f.size; i++)
    {
        f1 += f[i] * conj(f[i]);
    }
    return f1.real();
};
complex<double> vector_p(lattice_fermi r1, lattice_fermi r2)
{
    complex<double> ro(0.0, 0.0);
    for (int i = 0; i < r1.size; i++)
    {
        ro += (conj(r1[i]) * r2[i]);
    }
    return ro;
};
void Dslash2(lattice_fermi src, lattice_fermi &dest, lattice_gauge U, const double mass, const bool dag)
{
    dest.clean();
    const double a = 2.0;
    const complex<double> i(0.0, 1.0);
    complex<double> tmp;
    const double Half = 0.5;
    double flag = (dag == true) ? -1 : 1;
    for (int x = 0; x < src.lat_x; x++)
        for (int t = 0; t < src.lat_t; t++)
        {
            for (int s = 0; s < src.lat_spin; s++)
            {
                dest[(x * src.lat_t + t) * 2 + s] += -(a + mass) * src[(x * src.lat_t + t) * 2 + s];
            }
            int b_x = (x + src.lat_x - 1) % src.lat_x;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] + flag * src[(x * src.lat_t + t) * 2 + 1]) * Half * U[(b_x * src.lat_t + t) * 2 + 0];
            dest[(b_x * src.lat_t + t) * 2 + 0] += tmp;
            dest[(b_x * src.lat_t + t) * 2 + 1] += flag * tmp;
            int f_x = (x + 1) % src.lat_x;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] - flag * src[(x * src.lat_t + t) * 2 + 1]) * Half * conj(U[(x * src.lat_t + t) * 2 + 0]);
            dest[(f_x * src.lat_t + t) * 2 + 0] += tmp;
            dest[(f_x * src.lat_t + t) * 2 + 1] -= flag * tmp;
            int b_t = (t + src.lat_t - 1) % src.lat_t;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] + flag * i * src[(x * src.lat_t + t) * 2 + 1]) * Half * U[(x * src.lat_t + b_t) * 2 + 1];
            dest[(x * src.lat_t + b_t) * 2 + 0] += tmp;
            dest[(x * src.lat_t + b_t) * 2 + 1] -= flag * i * tmp;
            int f_t = (t + 1) % src.lat_t;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] - flag * i * src[(x * src.lat_t + t) * 2 + 1]) * Half * conj(U[(x * src.lat_t + t) * 2 + 1]);
            dest[(x * src.lat_t + f_t) * 2 + 0] += tmp;
            dest[(x * src.lat_t + f_t) * 2 + 1] += flag * i * tmp;
        }
}
void fermi_to_prop(lattice_fermi dest, lattice_propagator &prop, int i)
{
    for (int x = 0; x < dest.lat_x; x++)
        for (int t = 0; t < dest.lat_t; t++)
            for (int s = 0; s < dest.lat_spin; s++)
            {
                prop[(x * prop.lat_t + t) * prop.lat_spin * prop.lat_spin + (i * prop.lat_spin) + s] = dest[(x * prop.lat_t + t) * prop.lat_spin + s];
            }
}
int CG(lattice_fermi src, lattice_fermi &dest, lattice_gauge U, const double mass, const int max)
{
    lattice_fermi r0(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi rr0(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi z0(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi r1(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi z1(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi q(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi qq(src.lat_x, src.lat_t, src.lat_spin);
    lattice_fermi P(src.lat_x, src.lat_t, src.lat_spin);
    complex<double> aphi;
    complex<double> beta;
    for (int i = 0; i < src.size; i++)
    {
        r0[i] = 0;
    }
    for (int i = 0; i < dest.size; i++)
    {
        dest[i] = 0;
    }
    Dslash2(dest, rr0, U, mass, false);
    Dslash2(rr0, r0, U, mass, true);
    for (int i = 0; i < src.size; i++)
    {
        r0[i] = src[i] - r0[i];
    }
    for (int f = 1; f < max; f++)
    {
        std::complex<double> rho;
        rho = vector_p(r0, r0);
        std::complex<double> rho1;
        rho1 = vector_p(r1, r1);
        for (int i = 0; i < z0.size; i++)
        {
            z0[i] = r0[i];
        }
        if (f == 1)
        {
            for (int i = 0; i < P.size; i++)
            {
                P[i] = z0[i];
            }
        }
        else
        {
            beta = rho / rho1;
            for (int i = 0; i < P.size; i++)
                P[i] = z0[i] + beta * P[i];
        }
        Dslash2(P, qq, U, mass, false);
        Dslash2(qq, q, U, mass, true);
        aphi = rho / vector_p(P, q);
        for (int i = 0; i < dest.size; i++)
            dest[i] = dest[i] + aphi * P[i];
        for (int i = 0; i < r1.size; i++)
            r1[i] = r0[i];
        for (int i = 0; i < r0.size; i++)
            r0[i] = r0[i] - aphi * q[i];
        if (norm_2(r0) < 1e-12)
        {
            return 0;
        }
    }
    return 0;
}

class biCGstabX
// Stable Biconjugate Gradient Method
{
public:
    biCGstabX(int max4loop = 10,
              int nx = 4,
              int nt = 4,
              int ns = 2,
              float mass = 0.1,
              float min4diff = 1e-4,
              MatrixXcd *poit4A = NULL,
              VectorXcd *poit4b = NULL)
    {
        this->max4loop = max4loop;
        this->nu4dim = nx * nt * ns;
        this->nx = nx;
        this->nt = nt;
        this->ns = ns;
        this->mass = mass;
        this->min4diff = min4diff;
        this->poit4A = poit4A;
        this->poit4b = poit4b;
    };
    struct values2out
    {
        float diff;
        int loop;
        VectorXcd xi;
    };
    values2out _values2out;
    void Run()
    {
        Config();
        while (loop + 1)
        {
            Calculate();
            loop += 1;
            pi = pi1;
            xi = xi1;
            ri = ri1;
            proi = proi1;
            vi = vi1;
            wi = wi1;
        };
    };
    void del4biCGstabX(biCGstabX *poit)
    {
        poit = NULL;
        delete poit;
    };
    void del4values2out(values2out *poit)
    {
        poit = NULL;
        delete poit;
    };

protected:
    int max4loop;
    int nu4dim;
    int nx;
    int nt;
    int ns;
    float mass;
    float min4diff;
    void *poit4A;
    void *poit4b;
    VectorXcd b;
    VectorXcd r0;
    VectorXcd pi;
    VectorXcd pi1;
    VectorXcd xi;
    VectorXcd xi1;
    VectorXcd ri;
    VectorXcd ri1;
    VectorXcd vi;
    VectorXcd vi1;
    VectorXcd t;
    VectorXcd s;
    VectorXcd tmp4Vector;
    VectorXcd tmp4VectorX;
    complex<double> alpha;
    complex<double> beta;
    complex<double> proi;
    complex<double> proi1;
    complex<double> wi;
    complex<double> wi1;
    complex<double> tmp4Complex;
    complex<double> tmp4ComplexX;
    float diff;
    int loop;
    void Config()
    {
        loop = 0;
        proi = 1.0;
        alpha = 1.0;
        wi = 1.0;
        tmp4Complex = 0;
        tmp4ComplexX = 0;
        b = VectorXcd::Random(nu4dim);
        xi = VectorXcd::Random(nu4dim);
        vi1 = VectorXcd::Zero(nu4dim);
        t = VectorXcd::Zero(nu4dim);
        tmp4Vector = VectorXcd::Zero(nu4dim);
        tmp4VectorX = VectorXcd::Zero(nu4dim);
        Dslash(xi, tmp4Vector);
        ri = b - tmp4Vector;
        pi = VectorXcd::Zero(nu4dim);
        r0 = ri;
        vi = VectorXcd::Zero(nu4dim);
    };
    void Compare()
    {
        Dslash(xi, tmp4Vector);
        _values2out.diff = (tmp4Vector - b).norm() / b.norm();
        if (loop == max4loop || _values2out.diff < min4diff)
        {
            _values2out.loop = loop;
            _values2out.xi = xi;
            cout
                << "r(" << loop << "):" << endl
                << ri << "\n*************************************\n"
                << "x(" << loop << "):" << endl
                << xi << "\n*************************************\n"
                << "diff(" << loop << "):" << endl
                << _values2out.diff << "\n*************************************\n";
            loop = -2;
        };
    };
    void Calculate()
    // Refer to https://zh.m.wikipedia.org/wiki/稳定双共轭梯度法.
    {
        cout << "proi(" << loop << "):" << endl
             << proi << "\n*************************************\n"
             << "a(" << loop << "):" << endl
             << alpha << "\n*************************************\n"
             << "w(" << loop << "):" << endl
             << wi << "\n*************************************\n"
             << "vi(" << loop << "):" << endl
             << vi << "\n*************************************\n"
             << "pi(" << loop << "):" << endl
             << pi << "\n*************************************\n"
             << "r(" << loop << "):" << endl
             << ri << "\n*************************************\n"
             << "x(" << loop << "):" << endl
             << xi << "\n*************************************\n";
        proi1 = r0.transpose() * ri;
        beta = (proi1 / proi) * (alpha / wi);
        pi1 = ri + beta * (pi - wi * vi);
        Dslash(pi1, vi1);
        // vi1 = A * pi1;
        tmp4Complex = r0.transpose() * vi1;
        alpha = proi1 / tmp4Complex;
        s = ri - alpha * vi1;
        Dslash(s, t);
        // t = A * s;
        tmp4Complex = t.transpose() * s;
        tmp4ComplexX = t.transpose() * t;
        wi1 = tmp4Complex / tmp4ComplexX;
        xi1 = xi + alpha * pi1 + wi1 * s;
        Compare();
        ri1 = s - wi1 * t;
    };
    void Dslash(VectorXcd &value4src, VectorXcd &value4dest)
    {
        lattice_fermi src(nx, nt, ns);
        lattice_fermi dest(nx, nt, ns);
        lattice_gauge U(nx, nt, ns);
        for (int i = 0; i < U.size; i++)
        {
            U[i] = 1.0;
        };
        for (int i = 0; i < src.size; i++)
        {
            src[i] = value4src(i);
        };
        Dslash2(src, dest, U, mass, true);
        for (int i = 0; i < dest.size; i++)
        {
            value4dest(i) = dest[i];
        };
    };
};

int main()
{
    clock_t start = clock();


    biCGstabX bi;
    bi.Run();
    bi.del4biCGstabX(&bi);
    return 0;
}