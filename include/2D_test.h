#include <complex>
#include <vector>
#include <math.h>
using namespace std;
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