#include <iostream>
#include <complex>
#include <random>
#include <mpi.h>
#include <ctime>

class zvector
{
private:
public:
    int size;
    std::complex<double> *data;
    std::vector<int> dimensions;
    template <typename... Args>
    zvector(Args... args) : dimensions{args...}
    {
        size = 1;
        for (int d : dimensions)
        {
            size *= d;
        }
        data = new std::complex<double>[size];
    }
    zvector(const std::vector<int> &dimensions) : dimensions(dimensions)
    {
        size = 1;
        for (int d : dimensions)
        {
            size *= d;
        }
        data = new std::complex<double>[size];
    }

    ~zvector()
    {
        if (data != nullptr)
        {
            data = nullptr;
            delete[] data;
        }
    }

    template <typename... Args>
    std::complex<double> &operator()(Args... args)
    {
        int index = 0;
        int temp[] = {args...};
        int stride = 1;
        for (int i = sizeof...(args) - 1; i >= 0; --i)
        {
            index += temp[i] * stride;
            stride *= dimensions[i];
        }
        std::cout << "Position(1 dim): "
                  << index
                  << '\n';
        return data[index];
    }
    template <typename... Args>
    int dim1(Args... args)
    {
        int index = 0;
        int temp[] = {args...};
        int stride = 1;
        for (int i = sizeof...(args) - 1; i >= 0; --i)
        {
            index += temp[i] * stride;
            stride *= dimensions[i];
        }
        return index;
    }
    template <typename... Args>
    zvector *block(Args... args)
    {
        // int index = 0;
        std::vector<int> devide;
        devide = {args...};
        // int stride = 1;
        if (dimensions.size() != devide.size())
            throw std::invalid_argument("devide error.");
        int dsize = devide.size();
        int usize = size / dsize;
        zvector *list;
        std::complex<double> *tmp;
        std::cout << dsize << usize;
        // for (int i = 0; i < devide.size(); i++)
        // {
        //     for (int j = 0; j < usize; j++)
        //     {
        //        tmp[j] = data[i * usize + j];
        //     }
        //     // list[i].data = tmp;
        //     // list[i].display();
        // }
        return list;
    }

    zvector &operator=(const zvector &other)
    {
        if (this != &other)
        {
            if (size != other.size)
            {
                delete[] data;
                size = other.size;
                data = new std::complex<double>[size];
            }
            dimensions = other.dimensions;
            for (int i = 0; i < size; ++i)
            {
                data[i] = other.data[i];
            }
        }
        return *this;
    }

    template <typename T>
    zvector operator+(const T scalar) const
    {
        zvector result(dimensions);
        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] + scalar;
        }
        return result;
    }

    template <typename T>
    zvector operator-(const T scalar) const
    {
        zvector result(dimensions);
        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] - scalar;
        }
        return result;
    }

    template <typename T>
    zvector operator*(const T &scalar) const
    {
        zvector result(dimensions);
        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] * scalar;
        }
        return result;
    }

    template <typename T>
    zvector operator/(const T &scalar) const
    {
        if (scalar == 0)
            throw std::invalid_argument("Divisor cannot be zero.");

        zvector result(dimensions);
        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] / scalar;
        }

        return result;
    }
    zvector operator+(const zvector &other) const
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        zvector result(dimensions);

        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] + other.data[i];
        }

        return result;
    }

    zvector operator-(const zvector &other) const
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        zvector result(dimensions);

        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] - other.data[i];
        }

        return result;
    }

    zvector operator*(const zvector &other) const
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        zvector result(dimensions);

        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] * other.data[i];
        }

        return result;
    }

    zvector operator/(const zvector &other) const
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        for (int i = 0; i < other.size; ++i)
        {
            if (other.data[i] == std::complex<double>(0.0))
            {
                throw std::invalid_argument("Cannot divide by zero.");
            }
        }

        zvector result(dimensions);

        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] / other.data[i];
        }

        return result;
    }

    void assign_zero()
    {
        for (int i = 0; i < size; ++i)
            data[i] = 0;
    }

    void assign_unit()
    {
        for (int i = 0; i < size; ++i)
            data[i] = 1;
    }

    void assign_random()
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0.0, 1.0);
        for (int i = 0; i < size; ++i)
            data[i] = {dis(gen), dis(gen)};
    }

    double magnitude()
    {
        double sum = 0.0;
        for (int i = 0; i < size; ++i)
        {
            sum += std::norm(data[i]);
        }
        return sqrt(sum);
    }

    double norm2()
    {
        double sum = 0.0;
        for (int i = 0; i < size; ++i)
        {
            sum += std::norm(data[i]);
        }
        return sqrt(sum);
    }
    double norm2X()
    {
        double sum = 0.0;
        for (int i = 0; i < size; ++i)
        {
            sum += std::norm(data[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&sum, &sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return sqrt(sum);
    }

    zvector dot_product(const zvector &other)
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        zvector result(dimensions);

        for (int i = 0; i < size; ++i)
        {
            result.data[i] = data[i] * std::conj(other.data[i]);
        }

        return result;
    }

    std::complex<double> dot(const zvector &other)
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        std::complex<double> sum(0.0, 0.0);

        for (int i = 0; i < size; ++i)
        {
            sum += data[i] * std::conj(other.data[i]);
        }

        return sum;
    }
    std::complex<double> dotX(const zvector &other)
    {
        if (size != other.size)
        {
            throw std::invalid_argument("zvectors must have the same size.");
        }

        std::complex<double> sum(0.0, 0.0);

        for (int i = 0; i < size; ++i)
        {
            sum += data[i] * std::conj(other.data[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&sum, &sum, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);
        return sum;
    }
    zvector cross_product(const zvector &other)
    {
        if (size != 3 || other.size != 3)
        {
            throw std::invalid_argument("Cross product is only defined for zvectors of size 3.");
        }

        zvector result(dimensions);

        result.data[0] = data[1] * other.data[2] - data[2] * other.data[1];
        result.data[1] = data[2] * other.data[0] - data[0] * other.data[2];
        result.data[2] = data[0] * other.data[1] - data[1] * other.data[0];

        return result;
    }

    void display()
    {
        std::cout << "Dimensions: ";
        for (auto d : dimensions)
            std::cout << d << ' ';
        std::cout << '\n';

        std::cout << "Data: ";
        for (int i = 0; i < size && i < 10; ++i)
            std::cout << data[i] << ' ';

        if (size > 10)
        {
            std::cout << "...";
        }

        std::cout << '\n';
    }
};
class zdslash
{
public:
    zdslash(const zvector &U, int &lat_x, int &lat_t, int &lat_spin, double &mass, bool &dag) : U(U), lat_x(lat_x), lat_t(lat_t), lat_spin(lat_spin), mass(mass), dag(dag)
    {
    }
    ~zdslash()
    {
    }
    void _dslash(zvector &src, zvector &dest)
    {
        dest.assign_zero();
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
                    dest(x, t, s) += -(a + mass) * src(x, t, s);
                }

                // backward x
                int b_x = (x + lat_x - 1) % lat_x;
                tmp = (src(x, t, 0) + flag * src(x, t, 1)) * Half * U(b_x, t, 0);
                dest(b_x, t, 0) += tmp;
                dest(b_x, t, 1) += flag * tmp;

                // forward x
                int f_x = (x + 1) % lat_x;
                tmp = (src(x, t, 0) - flag * src(x, t, 1)) * Half * conj(U(x, t, 0));
                dest(f_x, t, 0) += tmp;
                dest(f_x, t, 1) -= flag * tmp;

                // backward t
                int b_t = (t + lat_t - 1) % lat_t;
                tmp = (src(x, t, 0) + flag * i * src(x, t, 1)) * Half * U(x, b_t, 1);
                dest(x, b_t, 0) += tmp;
                dest(x, b_t, 1) -= flag * i * tmp;

                // forward t
                int f_t = (t + 1) % lat_t;
                tmp = (src(x, t, 0) - flag * i * src(x, t, 1)) * Half * conj(U(x, t, 1));
                dest(x, f_t, 0) += tmp;
                dest(x, f_t, 1) += flag * i * tmp;
            }
    }

    void dslash(const zvector &src, zvector &dest)
    {
        dest = src * 0.2;
        dest = dest + 0.2;
    }

private:
    zvector U;
    int lat_x;
    int lat_t;
    int lat_spin;
    double mass;
    bool dag;
};
class zbicgstab
{
public:
    zbicgstab(int MAX_ITER, double &TOL, zvector &b, zdslash &Zdslash) : MAX_ITER(MAX_ITER), TOL(TOL), Zdslash(Zdslash), dimensions(b.dimensions)
    {
        this->b = &b;
    }

    ~zbicgstab()
    {
    }

    // solve the linear system Ax = b
    void solve()
    {
        zvector r = *b;//give b
        std::complex<double> rho_prev(1.0, 0.0);
        std::complex<double> rho(0.0, 0.0);
        std::complex<double> alpha(1.0, 0.0);
        std::complex<double> omega(1.0, 0.0);
        std::complex<double> beta(0.0, 0.0);
        zvector x(dimensions);
        zvector r_tilde(dimensions);
        zvector p(dimensions);
        zvector v(dimensions);
        zvector s(dimensions);
        zvector t(dimensions);

        x.assign_random(); // initial guess
        Zdslash.dslash(x, r_tilde);
        // zvector r = b - A * x;
        r = r - r_tilde;
        r_tilde = r;

        // x.assign_zero(); // if x=0;r_tilde=r0=b;
        // r_tilde = r;

        p.assign_zero();
        v.assign_zero();
        s.assign_zero();
        t.assign_zero();
        r.display();
        for (int i = 0; i < MAX_ITER; i++)
        {
            std::cout << "#"
                      << i
                      << "-Residual: "
                      << r.norm2X()
                      << std::endl;
            rho = r_tilde.dotX(r);
            beta = (rho / rho_prev) * (alpha / omega);
            p = r + (p - v * omega) * beta;
            // v = A * p;
            Zdslash.dslash(p, v);
            alpha = rho / r_tilde.dotX(v);
            s = r - v * alpha;
            // t = A * s;
            Zdslash.dslash(s, t);
            omega = t.dotX(s) / t.dotX(t);
            x = x + p * alpha + s * omega;
            r = s - t * omega;
            if (r.norm2X() < TOL)
            {
                std::cout << "#End-loop: "
                          << i
                          << std::endl;
                break;
            }
            rho_prev = rho;
        }
        std::cout << "#End-Residual: "
                  << r.norm2X()
                  << std::endl;
        x.display();
    }

private:
    int MAX_ITER;
    double TOL;
    zvector *b;
    zdslash Zdslash;
    std::vector<int> dimensions;
    // zvector block(zvector &src)
    // {
    //     zvector dest;
    //     int rank, size;
    //     MPI_Comm_size(MPI_COMM_WORLD, &size);
    //     MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //     ....
    // }
};
int main()
{
    MPI_Init(NULL, NULL);
    clock_t start = clock();
    int lat_x(100);
    int lat_t(100);
    int lat_spin(2); // const lat_s=2
    int MAX_ITER(1e6);
    double TOL(1e-6);
    zvector b(lat_x, lat_t, lat_spin);
    zvector U(lat_x, lat_t, lat_spin);
    b.assign_zero();
    b(0, 0, 0) = 1.0;
    b = b + 1e-18;
    U.assign_unit();
    double mass(1);
    bool dag(true);
    zdslash Zdslash(U, lat_x, lat_t, lat_spin, mass, dag);
    zbicgstab Zbicgstab(MAX_ITER, TOL, b, Zdslash);
    Zbicgstab.solve();
    clock_t end = clock();
    std::cout
        << "################"
        << "time cost:"
        << (double)(end - start) / CLOCKS_PER_SEC
        << "s"
        << std::endl;
    MPI_Finalize();
    return 0;
}
