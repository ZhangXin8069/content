#include "zx7.hpp"

BiCGstabX::BiCGstabX(int max4loop,
                     int dim4x,
                     int dim4t,
                     int dim4s,
                     double mass,
                     double min4diff)
{
    this->max4loop = max4loop;
    this->dim = dim4x * dim4t * dim4s;
    this->dim4x = dim4x;
    this->dim4t = dim4t;
    this->dim4s = dim4s;
    this->mass = mass;
    this->min4diff = min4diff;
};
BiCGstabX::~BiCGstabX(){
    // delete[] ri;
    // delete[] ri_prev;
    // delete[] xi;
    // delete[] xi_prev;
    // delete[] vi;
    // delete[] vi_prev;
    // delete[] tmp4vector;
    // delete[] tmp4vectorX;
};
void BiCGstabX::init2zero(complex<double> *&src1)
{











    for (int i = 0; i < dim; i++)
    {
        src1[i].real(0.0);
        src1[i].imag(0.0);
    };
};
void BiCGstabX::init2rand(complex<double> *&src1)
{
    for (int i = 0; i < dim; i++)
    {
        src1[i].real(rand() / 1e9);
        src1[i].imag(0.0);
    };
};
void BiCGstabX::init2one(complex<double> *&src1)
{

    for (int i = 0; i < dim; i++)
    {
        src1[i].real(1.0);
        src1[i].imag(0.0);
    };
};
void BiCGstabX::dot(complex<double> &dest, complex<double> *&src1, complex<double> *&src2)
{
    dest = (0.0, 0.0);
    for (int i = 0; i < dim; i++)
    {
        dest += src1[i] * src2[i];
    };
};
void BiCGstabX::mult(complex<double> *&dest, complex<double> &src1, complex<double> *&src2)
{
    for (int i = 0; i < dim; i++)
    {
        dest[i] = src1 * src2[i];
    };
};
void BiCGstabX::add(complex<double> *&dest, complex<double> *&src1, complex<double> *&src2)
{
    for (int i = 0; i < dim; i++)
    {
        dest[i] = src1[i] + src2[i];
    };
};
void BiCGstabX::cut(complex<double> *&dest, complex<double> *&src1, complex<double> *&src2)
{
    for (int i = 0; i < dim; i++)
    {
        dest[i] = src1[i] - src2[i];
    };
};
void BiCGstabX::equal(complex<double> *&src1, complex<double> *&src2)
{
    for (int i = 0; i < dim; i++)
    {
        src1[i] = src2[i];
    };
};
void BiCGstabX::print(complex<double> *&src1)
{
    for (int i = 0; i < dim; i++)
    {
        cout
            << "list(" << i << "):"
            << src1[i] << endl;
    };
};
void BiCGstabX::run()
{
    config();
    while (loop)
    {
        calculate();
        compare();
        reconfig();
    };
};
void BiCGstabX::reconfig()
{
    loop += 1;
    equal(pi_prev, pi);
    equal(xi_prev, xi);
    equal(ri_prev, ri);
    proi_prev = proi;
    equal(vi_prev, vi);
    wi_prev = wi;
};
void BiCGstabX::config()
{
    test = 2;
    loop = 1;
    proi = 1.0;
    proi_prev = 1.0;
    alpha = 1.0;
    wi = 1.0;
    wi_prev = 1.0;
    bool dag;
    dag = true;
    init2one(xi_prev);
    init2zero(vi_prev);
    init2zero(pi_prev);
    init2one(U);
    init2one(b);

    mult(tmp4vector, test, xi_prev);
    // dslash(xi_prev, tmp4vector);
    cut(r0, b, tmp4vector);
    equal(ri_prev, r0);
};
void BiCGstabX::compare()
{
    norm(diff, ri);
    cout
        << "diff=norm_r(" << loop << "):"
        << diff << endl;
    if (loop == max4loop || diff.real() < min4diff)
    {
        loop = -1;
    };
};
void BiCGstabX::calculate()
// Refer to https://zh.m.wikipedia.org/wiki/稳定双共轭梯度法.
{
    dot(proi, r0, ri_prev);
    tmp4complex = (proi / proi_prev) * (alpha / wi_prev);
    mult(tmp4vector, wi_prev, vi_prev);
    cut(tmp4vectorX, pi_prev, tmp4vector);
    mult(tmp4vector, tmp4complex, tmp4vectorX);
    add(pi, ri_prev, tmp4vector);

    // dslash(pi, vi);
    mult(vi, test, pi);

    dot(tmp4complex, r0, vi);
    alpha = proi / tmp4complex;
    mult(tmp4vector, alpha, vi);
    cut(s, ri_prev, tmp4vector);

    // dslash(s, tmp4vector);
    mult(tmp4vector, test, s);

    dot(tmp4complex, tmp4vector, s);
    dot(tmp4complexX, tmp4vector, tmp4vector);
    wi = tmp4complex / tmp4complexX;
    mult(tmp4vector, alpha, pi);
    mult(tmp4vectorX, wi, s);
    add(tmp4vector, tmp4vector, tmp4vectorX);
    add(xi, tmp4vector, xi_prev);
    mult(tmp4vector, wi, pi);
    cut(ri, s, tmp4vector);
    // print(ri);
    print(xi);
    // print(b);
};

void BiCGstabX::dslash(complex<double> *&src, complex<double> *&dest)
{
    const double a = 2.0;
    const complex<double> i(0.0, 1.0);
    const double Half = 0.5;
    double flag = (dag == true) ? -1 : 1;
    for (int src1 = 0; src1 < dim4x; src1++)
        for (int t = 0; t < dim4t; t++)
        {
            for (int s = 0; s < dim4s; s++)
            {
                dest[(src1 * dim4t + t) * 2 + s] += -(a + mass) * src[(src1 * dim4t + t) * 2 + s];
            };
            int b_x = (src1 + dim4x - 1) % dim4x;
            tmp4complex = (src[(src1 * dim4t + t) * 2 + 0] + flag * src[(src1 * dim4t + t) * 2 + 1]) * Half * U[(b_x * dim4t + t) * 2 + 0];
            dest[(b_x * dim4t + t) * 2 + 0] += tmp4complex;
            dest[(b_x * dim4t + t) * 2 + 1] += flag * tmp4complex;
            int f_x = (src1 + 1) % dim4x;
            tmp4complex = (src[(src1 * dim4t + t) * 2 + 0] - flag * src[(src1 * dim4t + t) * 2 + 1]) * Half * conj(U[(src1 * dim4t + t) * 2 + 0]);
            dest[(f_x * dim4t + t) * 2 + 0] += tmp4complex;
            dest[(f_x * dim4t + t) * 2 + 1] -= flag * tmp4complex;
            int b_t = (t + dim4t - 1) % dim4t;
            tmp4complex = (src[(src1 * dim4t + t) * 2 + 0] + flag * i * src[(src1 * dim4t + t) * 2 + 1]) * Half * U[(src1 * dim4t + b_t) * 2 + 1];
            dest[(src1 * dim4t + b_t) * 2 + 0] += tmp4complex;
            dest[(src1 * dim4t + b_t) * 2 + 1] -= flag * i * tmp4complex;
            int f_t = (t + 1) % dim4t;
            tmp4complex = (src[(src1 * dim4t + t) * 2 + 0] - flag * i * src[(src1 * dim4t + t) * 2 + 1]) * Half * conj(U[(src1 * dim4t + t) * 2 + 1]);
            dest[(src1 * dim4t + f_t) * 2 + 0] += tmp4complex;
            dest[(src1 * dim4t + f_t) * 2 + 1] += flag * i * tmp4complex;
        };
};

void BiCGstabX::norm(complex<double> &dest, complex<double> *&src1)
{
    dest = (0.0, 0.0);
    for (int i = 0; i < dim; i++)
    {
        dest += src1[i] * conj(src1[i]);
    };
};
int main()
{
    clock_t start = clock();
    BiCGstabX bi;
    bi.run();
    clock_t end = clock();
    cout
        // << "rank:"
        // << rank
        << "################"
        << "time cost:"
        << (double)(end - start) / CLOCKS_PER_SEC
        << "s"
        << endl;
    // MPI_Init(NULL, NULL);

    // MPI_Finalize();
}