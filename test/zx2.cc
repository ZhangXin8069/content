#include <iostream>
#include <complex>
#include "Eigen/Core"
using namespace std;
using namespace Eigen;

class biCGstab
// Stable Biconjugate Gradient Method
{
public:
    biCGstab(int max4loop = 10,
             int nu4dim = 2,
             float min4diff = 1e-5,
             MatrixXcd *poit4A = NULL,
             VectorXcd *poit4b = NULL)
    {
        this->max4loop = max4loop;
        this->nu4dim = nu4dim;
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
    void del4biCGstab(biCGstab *poit)
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
    float min4diff;
    void *poit4A;
    void *poit4b;
    MatrixXcd A;
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
        A = MatrixXcd::Random(nu4dim, nu4dim);
        b = VectorXcd::Random(nu4dim);
        xi = VectorXcd::Random(nu4dim);
        ri = b - A * xi;
        pi = VectorXcd::Zero(nu4dim);
        r0 = ri;
        vi = VectorXcd::Zero(nu4dim);
        cout << "A:" << endl
             << A
             << "\n*************************************\n"
             << "b:" << endl
             << b
             << "\n*************************************\n"
             << "x(0):"
             << endl
             << xi << endl
             << "\n*************************************\n"
             << "r(0):"
             << endl
             << ri << endl
             << "\n*************************************\n";
    };
    void Compare()
    {
        _values2out.diff = (A * xi - b).norm();
        cout
            << "diff(" << loop << "):" << endl
            << _values2out.diff << "\n*************************************\n";
        if (loop == max4loop || _values2out.diff < min4diff)
        {
            _values2out.loop = loop;
            _values2out.xi = xi;
            cout
                << "R(" << loop << "):" << endl
                << ri << "\n*************************************\n"
                << "X(" << loop << "):" << endl
                << xi << "\n*************************************\n";
            loop = -2;
        }
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
        vi1 = A * pi1;
        tmp4Complex = r0.transpose() * vi1;
        alpha = proi1 / tmp4Complex;
        s = ri - alpha * vi1;
        t = A * s;
        tmp4Complex = t.transpose() * s;
        tmp4ComplexX = t.transpose() * t;
        wi1 = tmp4Complex / tmp4ComplexX;
        xi1 = xi + alpha * pi1 + wi1 * s;
        Compare();
        ri1 = s - wi1 * t;
    };
};
int main()
{
    biCGstab bi(1e5, 1, 1e-5);
    bi.Run();
    biCGstab::values2out values = bi._values2out;
    cout << values.diff << endl;
    bi.del4values2out(&values);
    bi.del4biCGstab(&bi);
    return 0;
};
