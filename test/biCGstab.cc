// #include "biCGstab.h"
#include <iostream>
#include "Eigen/Core"
using namespace std;
using namespace Eigen;
class biCGstab
// Stable Biconjugate Gradient Method
{
public:
    biCGstab(int max4loop = 1e2,
             int nu4dim = 44,
             float min4diff = 1e-4,
             MatrixXd *poit4A = NULL,
             VectorXd *poit4b = NULL)
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
        VectorXd xi;
    };
    values2out _values2out;
    void Run()
    {
        Config();
        while (loop + 1)
        {
            Calculate();
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
    MatrixXd A;
    VectorXd b;
    VectorXd r0;
    VectorXd pi;
    VectorXd xi;
    VectorXd ri;
    VectorXd vi;
    double alpha;
    double proi;
    double wi;
    double tmp4double;
    float diff;
    int loop;
    void Config()
    {
        loop = 0;
        proi = 1.0;
        alpha = 1.0;
        wi = 1.0;
        A = MatrixXd::Random(nu4dim, nu4dim);
        b = VectorXd::Random(nu4dim);
        xi = VectorXd::Random(nu4dim);
        ri = b - A * xi;
        pi = VectorXd::Zero(nu4dim);
        r0 = ri;
        vi = VectorXd::Zero(nu4dim);
        cout << "A:" << endl
             << A << "\n*************************************\n"
             << "b:" << endl
             << b
             << "\n*************************************\n"
            //  << "x(0):"
            //  << endl
            //  << xi << endl
            //  << "\n*************************************\n"
            //  << "r(0):"
            //  << endl
            //  << ri << endl
            //  << "\n*************************************\n"
            ;
    };
    void Compare()
    {
        if (loop == max4loop || (A * xi - b).norm() < min4diff)
        {
            cout
                << "R(" << loop << "):" << endl
                << ri << "\n*************************************\n"
                << "X(" << loop << "):" << endl
                << xi << "\n*************************************\n";
            _values2out.loop = loop;
            _values2out.xi = xi;
            _values2out.diff = (A * xi - b).norm();
            loop = -1;
        }
    };
    void Calculate()
    // Refer to https://zh.m.wikipedia.org/wiki/稳定双共轭梯度法.
    {
        loop += 1;
        double &proI = tmp4double;
        proI = proi;

        double &proi1 = proi;
        proi1 = r0.adjoint() * ri;

        double &beta = tmp4double;
        beta = (proi1 / proI) * (alpha / wi);

        VectorXd &pi1 = pi;
        pi1 = ri + beta * (pi - wi * vi);

        VectorXd &vi1 = vi;
        vi1 = A * pi1; // dslash,A is special.

        alpha = proi1 / (r0.adjoint() * vi1);
        VectorXd &s = ri;
        s = ri - alpha * vi1;

        double &wi1 = wi;
        tmp4double = (A * s).adjoint() * s;
        wi1 = tmp4double / ((A * s).adjoint() * (A * s));

        VectorXd &xi1 = xi;
        xi1 = xi + alpha * pi1 + wi1 * s;
        Compare();
        VectorXd &ri1 = ri;
        ri1 = s - wi1 * (A * s);
        // cout << "w(" << loop << "):" << endl
        //      << wi << "\n*************************************\n"
        //      << "a(" << loop << "):" << endl
        //      << alpha << "\n*************************************\n"
        //      << "r(" << loop << "):" << endl
        //      << ri << "\n*************************************\n"
        //      << "x(" << loop << "):" << endl
        //      << xi << "\n*************************************\n";
    };
};
class biCGstabX : protected biCGstab
// Going to be Preconditioned Stable Biconjugate Gradient Method...
{
public:
    using biCGstab::_values2out;
    using biCGstab::biCGstab;
    using biCGstab::del4values2out;
    void del4biCGstab(biCGstabX *poit)
    {
        poit = NULL;
        delete poit;
    };
    void Run()
    {
        Config();
        while (loop + 1)
        {
            Calculate();
        };
    };

private:
    using biCGstab::Calculate;
    using biCGstab::Compare;
    using biCGstab::Config;
    void Calculate()
    // Refer to https://zh.m.wikipedia.org/wiki/稳定双共轭梯度法
    {
        loop += 1;
        double &proI = tmp4double;
        proI = proi;

        double &proi1 = proi;
        proi1 = r0.adjoint() * ri;

        double &beta = tmp4double;
        beta = (proi1 / proI) * (alpha / wi);

        VectorXd &pi1 = pi;
        pi1 = ri + beta * (pi - wi * vi);

        VectorXd &vi1 = vi;
        vi1 = A * pi1;

        alpha = proi1 / (r0.adjoint() * vi1);
        VectorXd &s = ri;
        s = ri - alpha * vi1;

        double &wi1 = wi;
        tmp4double = (A * s).adjoint() * s;
        wi1 = tmp4double / ((A * s).adjoint() * (A * s));

        VectorXd &xi1 = xi;
        xi1 = xi + alpha * pi1 + wi1 * s;
        Compare();
        VectorXd &ri1 = ri;
        ri1 = s - wi1 * (A * s);
        //     cout << "w(" << loop << "):" << endl
        //          << wi << "\n*************************************\n"
        //          << "a(" << loop << "):" << endl
        //          << alpha << "\n*************************************\n"
        //          << "r(" << loop << "):" << endl
        //          << ri << "\n*************************************\n"
        //          << "x(" << loop << "):" << endl
        //          << xi << "\n*************************************\n";
        //
    };
};