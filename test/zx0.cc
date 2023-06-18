#include "zx_h.h"
#include "mpi.h"
class biCGstabX
// Stable Biconjugate Gradient Method
{
public:
    biCGstabX(int max4loop = 10,
              int nx = 3,
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
    ~biCGstabX()
    {
        del4biCGstabX(this);
    }

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
    void del4biCGstabX(biCGstabX *poit)
    {
        poit = NULL;
        delete poit;
    };
};
int main()
{
    //     biCGstab bi;
    //     bi.Run();
    biCGstabX bi;
    bi.Run();
    return 0;
};