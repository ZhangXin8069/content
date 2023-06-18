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

            // mass term
            for (int s = 0; s < src.lat_spin; s++)
            {
                dest[(x * src.lat_t + t) * 2 + s] += -(a + mass) * src[(x * src.lat_t + t) * 2 + s];
            }

            // backward x
            int b_x = (x + src.lat_x - 1) % src.lat_x;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] + flag * src[(x * src.lat_t + t) * 2 + 1]) * Half * U[(b_x * src.lat_t + t) * 2 + 0];
            dest[(b_x * src.lat_t + t) * 2 + 0] += tmp;
            dest[(b_x * src.lat_t + t) * 2 + 1] += flag * tmp;

            // forward x
            int f_x = (x + 1) % src.lat_x;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] - flag * src[(x * src.lat_t + t) * 2 + 1]) * Half * conj(U[(x * src.lat_t + t) * 2 + 0]);
            dest[(f_x * src.lat_t + t) * 2 + 0] += tmp;
            dest[(f_x * src.lat_t + t) * 2 + 1] -= flag * tmp;

            // backward t
            int b_t = (t + src.lat_t - 1) % src.lat_t;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] + flag * i * src[(x * src.lat_t + t) * 2 + 1]) * Half * U[(x * src.lat_t + b_t) * 2 + 1];
            dest[(x * src.lat_t + b_t) * 2 + 0] += tmp;
            dest[(x * src.lat_t + b_t) * 2 + 1] -= flag * i * tmp;

            // forward t
            int f_t = (t + 1) % src.lat_t;
            tmp = (src[(x * src.lat_t + t) * 2 + 0] - flag * i * src[(x * src.lat_t + t) * 2 + 1]) * Half * conj(U[(x * src.lat_t + t) * 2 + 1]);
            dest[(x * src.lat_t + f_t) * 2 + 0] += tmp;
            dest[(x * src.lat_t + f_t) * 2 + 1] += flag * i * tmp;
        }
}