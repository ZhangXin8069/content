#include <complex>
#include <mpi.h>
#include <ctime>
#include <iostream>
using namespace std;
class BiCGstabX
// Stable Biconjugate Gradient Method
{
public:
    BiCGstabX(int max4loop = 100,
              int dim4x = 3,
              int dim4t = 3,
              int dim4s = 2,
              double mass = 0.1,
              double min4diff = 1e-5);
    ~BiCGstabX();
    void run();
    void config();
    void reconfig();
    void compare();
    void calculate();

private:
    int max4loop;
    int dim;
    int dim4x;
    int dim4t;
    int dim4s;
    double mass;
    double min4diff;
    int loop;
    complex<double> diff;
    bool dag;
    complex<double> test;
    complex<double> *U = new complex<double>[dim];
    complex<double> *pi = new complex<double>[dim];
    complex<double> *b = new complex<double>[dim];
    complex<double> *r0 = new complex<double>[dim];
    complex<double> *s = new complex<double>[dim];
    complex<double> *pi_prev = new complex<double>[dim];
    complex<double> *xi = new complex<double>[dim];
    complex<double> *xi_prev = new complex<double>[dim];
    complex<double> *ri = new complex<double>[dim];
    complex<double> *ri_prev = new complex<double>[dim];
    complex<double> *vi = new complex<double>[dim];
    complex<double> *vi_prev = new complex<double>[dim];
    complex<double> *tmp4vector = new complex<double>[dim];
    complex<double> *tmp4vectorX = new complex<double>[dim];
    complex<double> proi;
    complex<double> proi_prev;
    complex<double> wi;
    complex<double> alpha;
    complex<double> wi_prev;
    complex<double> tmp4complex;
    complex<double> tmp4complexX;

    void init2zero(complex<double> *&);
    void init2rand(complex<double> *&);
    void init2one(complex<double> *&);
    void dslash(complex<double> *&, complex<double> *&);
    void dot(complex<double> &, complex<double> *&, complex<double> *&);
    void mult(complex<double> *&, complex<double> &, complex<double> *&);
    void add(complex<double> *&, complex<double> *&, complex<double> *&);
    void cut(complex<double> *&, complex<double> *&, complex<double> *&);
    void equal(complex<double> *&, complex<double> *&);
    void print(complex<double> *&);
    void norm(complex<double> &, complex<double> *&);
};