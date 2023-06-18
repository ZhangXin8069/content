#ifndef _LATTICE_DSLASH_H
#define _LATTICE_DSLASH_H
#include "lattice_vec.hpp"
class lattice_dslash
{
public:
    lattice_dslash(lattice_vec &U, int &lat_x, int &lat_t, int &lat_spin, double &mass, bool &dag) : U(U), lat_x(lat_x), lat_t(lat_t), lat_spin(lat_spin), mass(mass), dag(dag)
    {
    }

    ~lattice_dslash()
    {
    }
    void dslash(lattice_vec &src, lattice_vec &dest)
    {
        dest.clean();
        double a = 2.0;
        std::complex<double> i(0.0, 1.0);
        std::complex<double> tmp;
        double Half = 0.5;
        double flag = (dag == true) ? -1 : 1;
        for (int x = 0; x < lat_x; x++)
            for (int t = 0; t < lat_t; t++)
            {
                // mass term
                for (int s = 0; s < lat_spin; s++)
                {
                    dest[(x * lat_t + t) * 2 + s] += -(a + mass) * src[(x * lat_t + t) * 2 + s];
                }

                // backward x
                int b_x = (x + lat_x - 1) % lat_x;
                tmp = (src[(x * lat_t + t) * 2 + 0] + flag * src[(x * lat_t + t) * 2 + 1]) * Half * U[(b_x * lat_t + t) * 2 + 0];
                dest[(b_x * lat_t + t) * 2 + 0] += tmp;
                dest[(b_x * lat_t + t) * 2 + 1] += flag * tmp;

                // forward x
                int f_x = (x + 1) % lat_x;
                tmp = (src[(x * lat_t + t) * 2 + 0] - flag * src[(x * lat_t + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t) * 2 + 0]);
                dest[(f_x * lat_t + t) * 2 + 0] += tmp;
                dest[(f_x * lat_t + t) * 2 + 1] -= flag * tmp;

                // backward t
                int b_t = (t + lat_t - 1) % lat_t;
                tmp = (src[(x * lat_t + t) * 2 + 0] + flag * i * src[(x * lat_t + t) * 2 + 1]) * Half * U[(x * lat_t + b_t) * 2 + 1];
                dest[(x * lat_t + b_t) * 2 + 0] += tmp;
                dest[(x * lat_t + b_t) * 2 + 1] -= flag * i * tmp;

                // forward t
                int f_t = (t + 1) % lat_t;
                tmp = (src[(x * lat_t + t) * 2 + 0] - flag * i * src[(x * lat_t + t) * 2 + 1]) * Half * conj(U[(x * lat_t + t) * 2 + 1]);
                dest[(x * lat_t + f_t) * 2 + 0] += tmp;
                dest[(x * lat_t + f_t) * 2 + 1] += flag * i * tmp;
            }
    }

    void _dslash(lattice_vec &src, lattice_vec &dest)
    {
        dest = src * 2.0;
        dest[0] = 0.5;
    }

private:
    lattice_vec U;
    int lat_x;
    int lat_t;
    int lat_spin;
    double mass;
    bool dag;
};

#endif