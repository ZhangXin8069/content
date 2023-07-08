#include <iostream>
#include <random>
#include <mpi.h>
using namespace std;
class Complex
{
public:
    double data[2];
    Complex(double r = 0.0, double i = 0.0)
    {
        data[0] = r;
        data[1] = i;
    }
    Complex operator+(const Complex &other) const
    {
        return Complex(data[0] + other.data[0], data[1] + other.data[1]);
    }
    Complex operator-(const Complex &other) const
    {
        return Complex(data[0] - other.data[0], data[1] - other.data[1]);
    }
    Complex operator*(const Complex &other) const
    {
        return Complex(data[0] * other.data[0] - data[1] * other.data[1],
                       data[0] * other.data[1] + data[1] * other.data[0]);
    }
    Complex operator*(const double &other) const
    {
        return Complex(data[0] * other, data[1] * other);
    }
    Complex operator/(const Complex &other) const
    {
        double denom = other.data[0] * other.data[0] + other.data[1] * other.data[1];
        return Complex((data[0] * other.data[0] + data[1] * other.data[1]) / denom,
                       (data[1] * other.data[0] - data[0] * other.data[1]) / denom);
    }
    Complex operator/(const double &other) const
    {
        return Complex(data[0] / other, data[1] / other);
    }
    Complex operator-() const
    {
        return Complex(-data[0], -data[1]);
    }
    Complex &operator+=(const Complex &other)
    {
        data[0] += other.data[0];
        data[1] += other.data[1];
        return *this;
    }
    Complex &operator-=(const Complex &other)
    {
        data[0] -= other.data[0];
        data[1] -= other.data[1];
        return *this;
    }
    Complex &operator*=(const Complex &other)
    {
        double real = data[0] * other.data[0] - data[1] * other.data[1];
        double imag = data[0] * other.data[1] + data[1] * other.data[0];
        data[0] = real;
        data[1] = imag;
        return *this;
    }
    Complex &operator*=(const double &scalar)
    {
        data[0] *= scalar;
        data[1] *= scalar;
        return *this;
    }
    Complex &operator/=(const Complex &other)
    {
        double denom = other.data[0] * other.data[0] + other.data[1] * other.data[1];
        double real = (data[0] * other.data[0] + data[1] * other.data[1]) / denom;
        double imag = (data[1] * other.data[0] - data[0] * other.data[1]) / denom;
        data[0] = real;
        data[1] = imag;
        return *this;
    }
    Complex &operator/=(const double &other)
    {
        data[0] /= other;
        data[1] /= other;
        return *this;
    }
    bool operator==(const Complex &other) const
    {
        return (data[0] == other.data[0] && data[1] == other.data[1]);
    }
    bool operator!=(const Complex &other) const
    {
        return !(*this == other);
    }
    friend std::ostream &operator<<(std::ostream &os, const Complex &c)
    {
        if (c.data[1] >= 0.0)
        {
            os << c.data[0] << " + " << c.data[1] << "i";
        }
        else
        {
            os << c.data[0] << " - " << std::abs(c.data[1]) << "i";
        }
        return os;
    }
    Complex conj()
    {
        return Complex(data[0], data[1]);
    }
};

class LatticeFermi
{
public:
    int lat_x, lat_y, lat_z, lat_t, lat_s, lat_c;
    int size;
    Complex *lattice_vec;
    LatticeFermi(int lat_x, int lat_y, int lat_z, int lat_t, int lat_s, int lat_c)
        : lat_x(lat_x), lat_y(lat_y), lat_z(lat_z), lat_t(lat_t), lat_s(lat_s), lat_c(lat_c), size(lat_x * lat_y * lat_z * lat_t * lat_s * lat_c)
    {
        lattice_vec = new Complex[size];
    }
    ~LatticeFermi()
    {
        if (lattice_vec != nullptr)
        {
            lattice_vec = nullptr;
            delete[] lattice_vec;
        }
    }
    void assign_zero()
    {
        for (int i = 0; i < size; i++)
        {
            lattice_vec[i].data[0] = 0;
            lattice_vec[i].data[1] = 0;
        }
    }
    void assign_unit()
    {
        for (int i = 0; i < size; i++)
        {
            lattice_vec[i].data[0] = 1;
            lattice_vec[i].data[1] = 0;
        }
    }
    void assign_random(unsigned seed = 32767)
    {
        std::default_random_engine e(seed);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < size; i++)
        {
            lattice_vec[i].data[0] = u(e);
            lattice_vec[i].data[1] = u(e);
        }
    }
    void info()
    {
        std::cout << "lat_x:" << lat_x << std::endl;
        std::cout << "lat_y:" << lat_y << std::endl;
        std::cout << "lat_z:" << lat_z << std::endl;
        std::cout << "lat_t:" << lat_t << std::endl;
        std::cout << "lat_s:" << lat_s << std::endl;
        std::cout << "lat_c:" << lat_c << std::endl;
        std::cout << "size:" << size << std::endl;
    }
    const Complex &operator[](int index) const
    {
        return lattice_vec[index];
    }
    Complex &operator[](int index)
    {
        return lattice_vec[index];
    }
    const Complex &operator()(int index_x, int index_y, int index_z, int index_t, int index_s, int index_c) const
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c + index_y * lat_z * lat_t * lat_s * lat_c + index_z * lat_t * lat_s * lat_c + index_t * lat_s * lat_c + index_s * lat_c + index_c;
        return lattice_vec[index];
    }
    Complex &operator()(int index_x, int index_y, int index_z, int index_t, int index_s, int index_c)
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c + index_y * lat_z * lat_t * lat_s * lat_c + index_z * lat_t * lat_s * lat_c + index_t * lat_s * lat_c + index_s * lat_c + index_c;
        return lattice_vec[index];
    }
    LatticeFermi operator+(const LatticeFermi &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] + other[i];
        }
        return result;
    }
    LatticeFermi operator-(const LatticeFermi &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] - other[i];
        }
        return result;
    }
    LatticeFermi operator-() const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result[i] = -lattice_vec[i];
        }
        return result;
    }
    LatticeFermi operator*(const LatticeFermi &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] * other[i];
        }
        return result;
    }
    LatticeFermi operator/(const LatticeFermi &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] / other[i];
        }
        return result;
    }
    LatticeFermi operator+(const Complex &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] + other;
        }
        return result;
    }
    LatticeFermi operator-(const Complex &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] - other;
        }
        return result;
    }
    LatticeFermi operator*(const Complex &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] * other;
        }
        return result;
    }
    LatticeFermi operator/(const Complex &other) const
    {
        LatticeFermi result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] / other;
        }
        return result;
    }
    bool operator==(const LatticeFermi &other) const
    {
        if (size != other.size)
        {
            return false;
        }
        for (int i = 0; i < size; ++i)
        {
            if (lattice_vec[i] != other[i])
            {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const LatticeFermi &other) const
    {
        return !(*this == other);
    }
    void print(int index_x, int index_y, int index_z, int index_t, int index_s, int index_c)
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c + index_y * lat_z * lat_t * lat_s * lat_c + index_z * lat_t * lat_s * lat_c + index_t * lat_s * lat_c + index_s * lat_c + index_c;
        std::cout << "lattice_vec[" << index_x << "][" << index_y << "][" << index_z << "][" << index_t << "][" << index_s << "][" << index_c << "] = " << lattice_vec[index] << std::endl;
    }
    void print()
    {
        for (int x = 0; x < lat_x; x++)
        {
            for (int y = 0; y < lat_y; y++)
            {
                for (int z = 0; z < lat_z; z++)
                {
                    for (int t = 0; t < lat_t; t++)
                    {
                        for (int s = 0; s < lat_s; s++)
                        {
                            for (int c = 0; c < lat_c; c++)
                            {
                                print(x, y, z, t, s, c);
                            }
                        }
                    }
                }
            }
        }
    }
    double norm_2()
    {
        double result = 0;
        for (int i = 0; i < size; i++)
        {
            result = result + lattice_vec[i].data[0] * lattice_vec[i].data[0] + lattice_vec[i].data[1] * lattice_vec[i].data[1];
        }
        return result;
    }
    double norm_2X()
    {
        double local_result = 0;
        double global_result = 0;
        local_result = norm_2();
        MPI_Allreduce(&local_result, &global_result, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return global_result;
    }
    Complex dot(const LatticeFermi &other)
    {
        Complex result;
        for (int i = 0; i < size; i++)
        {
            result = result + lattice_vec[i].conj() * other[i];
        }
        return result;
    }
    Complex dotX(const LatticeFermi &other)
    {
        Complex local_result;
        Complex global_result;
        local_result = dot(other);
        MPI_Allreduce(&local_result, &global_result, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return global_result;
    }
    LatticeFermi block(int num_x, int num_y, int num_z, int num_t, int index_x, int index_y, int index_z, int index_t)
    {
        int block_x, block_y, block_z, block_t;
        block_x = lat_x / num_x;
        block_y = lat_y / num_y;
        block_z = lat_z / num_z;
        block_t = lat_t / num_t;
        int start_x = index_x * block_x;
        int start_y = index_y * block_y;
        int start_z = index_z * block_z;
        int start_t = index_t * block_t;
        LatticeFermi result(block_x, block_y, block_z, block_t, lat_s, lat_c);
        int global_x, global_y, global_z, global_t, global_index;
        for (int x = 0; x < block_x; x++)
        {
            global_x = start_x + x;
            for (int y = 0; y < block_y; y++)
            {
                global_y = start_y + y;
                for (int z = 0; z < block_z; z++)
                {
                    global_z = start_z + z;
                    for (int t = 0; t < block_t; t++)
                    {
                        global_t = start_t + t;
                        for (int s = 0; s < lat_s; s++)
                        {
                            for (int c = 0; c < lat_c; c++)
                            {
                                global_index = global_x * lat_y * lat_z * lat_t * lat_s * lat_c + global_y * lat_z * lat_t * lat_s * lat_c + global_z * lat_t * lat_s * lat_c + global_t * lat_s * lat_c + s * lat_c + c;
                                result(x, y, z, t, s, c) = lattice_vec[global_index];
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
    LatticeFermi block(int num_x, int num_y, int num_z, int num_t)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return block(num_x, num_y, num_z, num_t,
                     rank / (num_y * num_z * num_t),
                     (rank / (num_z * num_t)) % num_y,
                     (rank / num_t) % num_z,
                     rank % num_t);
    }
    LatticeFermi reback(int num_x, int num_y, int num_z, int num_t)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int start_x = (rank / (num_y * num_z * num_t)) * lat_x;
        int start_y = ((rank / (num_z * num_t)) % num_y) * lat_y;
        int start_z = ((rank / num_t) % num_z) * lat_z;
        int start_t = (rank % num_t) * lat_t;
        LatticeFermi result(num_x * lat_x, num_y * lat_y, num_z * lat_z, num_t * lat_t, lat_s, lat_c);
        int global_x, global_y, global_z, global_t, index;
        for (int x = 0; x < lat_x; x++)
        {
            global_x = start_x + x;
            for (int y = 0; y < lat_y; y++)
            {
                global_y = start_y + y;
                for (int z = 0; z < lat_z; z++)
                {
                    global_z = start_z + z;
                    for (int t = 0; t < lat_t; t++)
                    {
                        global_t = start_t + t;
                        for (int s = 0; s < lat_s; s++)
                        {
                            for (int c = 0; c < lat_c; c++)
                            {
                                index = x * lat_y * lat_z * lat_t * lat_s * lat_c + y * lat_z * lat_t * lat_s * lat_c + z * lat_t * lat_s * lat_c + t * lat_s * lat_c + s * lat_c + c;
                                result(global_x, global_y, global_z, global_t, s, c) = lattice_vec[index];
                            }
                        }
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, result.lattice_vec, result.size * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return result;
    }
};
class Gamme
{
    /*
    Gamme0=
    [[0,0,0,i],
    [0,0,i,0],
    [0,-i,0,0],
    [-i,0,0,0]]
    Gamme1=
    [[0,0,0,-1],
    [0,0,1,0],
    [0,1,0,0],
    [-1,0,0,0]]
    Gamme2=
    [[0,0,i,0],
    [0,0,0,-i],
    [-i,0,0,0],
    [0,i,0,0]]
    Gamme3=
    [[0,0,1,0],
    [0,0,0,1],
    [1,0,0,0],
    [0,1,0,0]]
    */
};

class LatticeGauge
{
public:
    int lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0, lat_c1;
    int size;
    Complex *lattice_vec;
    LatticeGauge(int lat_x, int lat_y, int lat_z, int lat_t, int lat_s, int lat_c)
        : lat_x(lat_x), lat_y(lat_y), lat_z(lat_z), lat_t(lat_t), lat_s(lat_s), lat_c0(lat_c), lat_c1(lat_c), size(lat_x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1)
    {
        lattice_vec = new Complex[size];
    }
    ~LatticeGauge()
    {
        if (lattice_vec != nullptr)
        {
            lattice_vec = nullptr;
            delete[] lattice_vec;
        }
    }
    void assign_zero()
    {
        for (int i = 0; i < size; i++)
        {
            lattice_vec[i].data[0] = 0;
            lattice_vec[i].data[1] = 0;
        }
    }
    void assign_unit()
    {
        for (int i = 0; i < size; i++)
        {
            lattice_vec[i].data[0] = 1;
            lattice_vec[i].data[1] = 0;
        }
    }
    void assign_random(unsigned seed = 32767)
    {
        std::default_random_engine e(seed);
        std::uniform_real_distribution<double> u(0.0, 1.0);
        for (int i = 0; i < size; i++)
        {
            lattice_vec[i].data[0] = u(e);
            lattice_vec[i].data[1] = u(e);
        }
    }
    void info()
    {
        std::cout << "lat_x:" << lat_x << std::endl;
        std::cout << "lat_y:" << lat_y << std::endl;
        std::cout << "lat_z:" << lat_z << std::endl;
        std::cout << "lat_t:" << lat_t << std::endl;
        std::cout << "lat_s:" << lat_s << std::endl;
        std::cout << "lat_c0:" << lat_c0 << std::endl;
        std::cout << "lat_c1:" << lat_c1 << std::endl;
        std::cout << "size:" << size << std::endl;
    }
    const Complex &operator[](int index) const
    {
        return lattice_vec[index];
    }
    Complex &operator[](int index)
    {
        return lattice_vec[index];
    }
    const Complex &operator()(int index_x, int index_y, int index_z, int index_t, int index_s, int index_c0, int index_c1) const
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_z * lat_t * lat_s * lat_c0 * lat_c1 + index_t * lat_s * lat_c0 * lat_c1 + index_s * lat_c0 * lat_c1 + index_c0 * lat_c1 + index_c1;
        return lattice_vec[index];
    }
    Complex &operator()(int index_x, int index_y, int index_z, int index_t, int index_s, int index_c0, int index_c1)
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_z * lat_t * lat_s * lat_c0 * lat_c1 + index_t * lat_s * lat_c0 * lat_c1 + index_s * lat_c0 * lat_c1 + index_c0 * lat_c1 + index_c1;
        return lattice_vec[index];
    }
    LatticeGauge operator+(const LatticeGauge &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] + other[i];
        }
        return result;
    }
    LatticeGauge operator-(const LatticeGauge &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] - other[i];
        }
        return result;
    }
    LatticeGauge operator-() const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result[i] = -lattice_vec[i];
        }
        return result;
    }
    LatticeGauge operator*(const LatticeGauge &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] * other[i];
        }
        return result;
    }
    LatticeGauge operator/(const LatticeGauge &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result[i] = lattice_vec[i] / other[i];
        }
        return result;
    }
    LatticeGauge operator+(const Complex &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] + other;
        }
        return result;
    }
    LatticeGauge operator-(const Complex &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] - other;
        }
        return result;
    }
    LatticeGauge operator*(const Complex &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] * other;
        }
        return result;
    }
    LatticeGauge operator/(const Complex &other) const
    {
        LatticeGauge result(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c0);
        for (int i = 0; i < size; ++i)
        {
            result.lattice_vec[i] = lattice_vec[i] / other;
        }
        return result;
    }
    bool operator==(const LatticeGauge &other) const
    {
        if (size != other.size)
        {
            return false;
        }
        for (int i = 0; i < size; ++i)
        {
            if (lattice_vec[i] != other[i])
            {
                return false;
            }
        }
        return true;
    }
    bool operator!=(const LatticeGauge &other) const
    {
        return !(*this == other);
    }
    void print(int index_x, int index_y, int index_z, int index_t, int index_s, int index_c0, int index_c1)
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_z * lat_t * lat_s * lat_c0 * lat_c1 + index_t * lat_s * lat_c0 * lat_c1 + index_s * lat_c0 * lat_c1 + index_c0 * lat_c1 + index_c1;
        std::cout << "lattice_vec[" << index_x << "][" << index_y << "][" << index_z << "][" << index_t << "][" << index_s << "][" << index_c0 << "][" << index_c1 << "] = " << lattice_vec[index] << std::endl;
    }
    void print()
    {
        for (int x = 0; x < lat_x; x++)
        {
            for (int y = 0; y < lat_y; y++)
            {
                for (int z = 0; z < lat_z; z++)
                {
                    for (int t = 0; t < lat_t; t++)
                    {
                        for (int s = 0; s < lat_s; s++)
                        {
                            for (int c0 = 0; c0 < lat_c0; c0++)
                            {
                                for (int c1 = 0; c1 < lat_c1; c1++)
                                {
                                    print(x, y, z, t, s, c0, c1);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    LatticeGauge block(int num_x, int num_y, int num_z, int num_t, int index_x, int index_y, int index_z, int index_t)
    {
        int block_x, block_y, block_z, block_t;
        block_x = lat_x / num_x;
        block_y = lat_y / num_y;
        block_z = lat_z / num_z;
        block_t = lat_t / num_t;
        int start_x = index_x * block_x;
        int start_y = index_y * block_y;
        int start_z = index_z * block_z;
        int start_t = index_t * block_t;
        LatticeGauge result(block_x, block_y, block_z, block_t, lat_s, lat_c0);
        int global_x, global_y, global_z, global_t, global_index;
        for (int x = 0; x < block_x; x++)
        {
            global_x = start_x + x;
            for (int y = 0; y < block_y; y++)
            {
                global_y = start_y + y;
                for (int z = 0; z < block_z; z++)
                {
                    global_z = start_z + z;
                    for (int t = 0; t < block_t; t++)
                    {
                        global_t = start_t + t;
                        for (int s = 0; s < lat_s; s++)
                        {
                            for (int c0 = 0; c0 < lat_c0; c0++)
                            {
                                for (int c1 = 0; c1 < lat_c1; c1++)
                                {
                                    global_index = global_x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + global_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + global_z * lat_t * lat_s * lat_c0 * lat_c1 + global_t * lat_s * lat_c0 * lat_c1 + s * lat_c0 * lat_c1 + c0 * lat_c1 + c1;
                                    result(x, y, z, t, s, c0, c1) = lattice_vec[global_index];
                                }
                            }
                        }
                    }
                }
            }
        }
        return result;
    }
    LatticeGauge block(int num_x, int num_y, int num_z, int num_t)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return block(num_x, num_y, num_z, num_t,
                     rank / (num_y * num_z * num_t),
                     (rank / (num_z * num_t)) % num_y,
                     (rank / num_t) % num_z,
                     rank % num_t);
    }
    LatticeGauge reback(int num_x, int num_y, int num_z, int num_t)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        int start_x = (rank / (num_y * num_z * num_t)) * lat_x;
        int start_y = ((rank / (num_z * num_t)) % num_y) * lat_y;
        int start_z = ((rank / num_t) % num_z) * lat_z;
        int start_t = (rank % num_t) * lat_t;
        LatticeGauge result(num_x * lat_x, num_y * lat_y, num_z * lat_z, num_t * lat_t, lat_s, lat_c0);
        int global_x, global_y, global_z, global_t, index;
        for (int x = 0; x < lat_x; x++)
        {
            global_x = start_x + x;
            for (int y = 0; y < lat_y; y++)
            {
                global_y = start_y + y;
                for (int z = 0; z < lat_z; z++)
                {
                    global_z = start_z + z;
                    for (int t = 0; t < lat_t; t++)
                    {
                        global_t = start_t + t;
                        for (int s = 0; s < lat_s; s++)
                        {
                            for (int c0 = 0; c0 < lat_c0; c0++)
                            {
                                for (int c1 = 0; c1 < lat_c1; c1++)
                                {
                                    index = x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + z * lat_t * lat_s * lat_c0 * lat_c1 + t * lat_s * lat_c0 * lat_c1 + s * lat_c0 * lat_c1 + c0 * lat_c1 + c1;
                                    result(global_x, global_y, global_z, global_t, s, c0, c1) = lattice_vec[index];
                                }
                            }
                        }
                    }
                }
            }
        }
        MPI_Allreduce(MPI_IN_PLACE, result.lattice_vec, result.size * 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        return result;
    }
};
void dslash(LatticeGauge &U, LatticeFermi &src, LatticeFermi &dest, int num_x, int num_y, int num_z, int num_t, double mass = 1.0, bool dag = true)
{
    dest.assign_zero();
    const double a = 2.0;
    const Complex i(0.0, 1.0);
    Complex tmp;
    Complex tmp0;
    Complex tmp1;
    const double Half = 0.5;
    double flag = (dag == true) ? -1 : 1;
    for (int x = 0; x < src.lat_x; x++)
    {
        for (int y = 0; y < src.lat_y; y++)
        {
            for (int z = 0; z < src.lat_z; z++)
            {
                for (int t = 0; t < src.lat_t; t++)
                {
                    int c = 0;
                    // mass term
                    for (int s = 0; s < src.lat_s; s++)
                    {
                        dest(x, y, z, t, s, c) += -src(x, y, z, t, s, c) * (a + mass);
                    }
                    int s = 0;
                    // backward x
                    int b_x = (x + src.lat_x - 1) % src.lat_x;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // forward x
                    int f_x = (x + 1) % src.lat_x;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // backward y
                    int b_y = (y + src.lat_y - 1) % src.lat_y;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // forward y
                    int f_y = (y + 1) % src.lat_y;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // backward z
                    int b_z = (z + src.lat_z - 1) % src.lat_z;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // forward z
                    int f_z = (z + 1) % src.lat_z;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // backward t
                    int b_t = (t + src.lat_t - 1) % src.lat_t;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                    // forward t
                    int f_t = (t + 1) % src.lat_t;
                    tmp = (src(x, y, z, t, s, c) - src(x, y, z, t, s, c)) * U(x, y, z, t, 1, 0, 0) * flag * Half;
                    dest(x, y, z, t, 0, 0) += tmp;
                    dest(x, y, z, t, 1, 0) += tmp * flag;
                    dest(x, y, z, t, 2, 0) += tmp;
                    dest(x, y, z, t, 3, 0) += tmp * flag;
                }
            }
        }
    }
}
void dslash(LatticeGauge &U, LatticeFermi &src, LatticeFermi &dest, double mass = 1.0, bool dag = true)
{
    dest = src * 0.5;
}
void cg(LatticeGauge &U0, LatticeFermi &b0, LatticeFermi &x0, int num_x, int num_y, int num_z, int num_t, double mass = 1.0, bool dag = true, int MAX_ITER = 1e6, double TOL = 1e-6)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Complex rho_prev(1.0, 0.0), rho(0.0, 0.0), alpha(1.0, 0.0), omega(1.0, 0.0), beta(0.0, 0.0);
    LatticeGauge U = U0.block(num_x, num_y, num_z, num_t);
    LatticeFermi b = b0.block(num_x, num_y, num_z, num_t);
    x0.assign_zero();
    LatticeFermi x = x0.block(num_x, num_y, num_z, num_t);
    LatticeFermi
        r = b,
        r_tilde = b,
        p(b.lat_x, b.lat_y, b.lat_z, b.lat_t, b.lat_s, b.lat_c),
        v(b.lat_x, b.lat_y, b.lat_z, b.lat_t, b.lat_s, b.lat_c),
        s(b.lat_x, b.lat_y, b.lat_z, b.lat_t, b.lat_s, b.lat_c),
        t(b.lat_x, b.lat_y, b.lat_z, b.lat_t, b.lat_s, b.lat_c);
    // x.rand(); // initial guess
    // dslash(x, r_tilde);
    // // ComplexVector r = b_ - A * x;
    // r = r - r_tilde;
    // if x=0;r_tilde = r0 = b_;
    for (int i = 0; i < MAX_ITER; i++)
    {
        rho = r_tilde.dotX(r);
        beta = (rho / rho_prev) * (alpha / omega);
        p = r + (p - v * omega) * beta;
        // v = A * p;
        // dslash(U, p, v, num_x, num_y, num_z, num_t, mass, dag);
        dslash(U, p, v);
        alpha = rho / r_tilde.dotX(v);
        s = r - v * alpha;
        // t = A * s;
        // dslash(U, s, t, num_x, num_y, num_z, num_t, mass, dag);
        dslash(U, s, t);
        omega = t.dotX(s) / t.dotX(t);
        x = x + p * alpha + s * omega;
        r = s - t * omega;
        if (r.norm_2X() < TOL || i == MAX_ITER - 1)
        {
            std::cout << "##loop "
                      << i
                      << "##Residual:"
                      << r.norm_2X()
                      << std::endl;
            break;
        }
        rho_prev = rho;
    }
}
int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    double start, end;
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    start = MPI_Wtime();
    int lat_x(16), lat_y(16), lat_z(16), lat_t(32), lat_s(4), lat_c(3);
    int num_x(1), num_y(2), num_z(4), num_t(2);
    int MAX_ITER(1e6);
    double TOL(1e-6);
    LatticeGauge U(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    LatticeFermi b(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    LatticeFermi x(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    U.assign_random(111);
    b.assign_random(222);
    x.assign_random(333);
    // cg(U, b, x, num_x, num_y, num_z, num_t);
    LatticeFermi block_b = b.block(num_x, num_y, num_z, num_t);
    LatticeFermi block_x = x.block(num_x, num_y, num_z, num_t);
    block_x.info();
    Complex dotxb = x.dot(b);
    Complex dotxb_ = block_x.dot(block_b);
    Complex dotXxb = block_x.dotX(block_b);
    cout << "dotxb_:" << dotxb_ << endl;
    cout << "dotxb:" << dotxb << endl;
    cout << "dotxb:" << dotxb << endl;
    cout << "dotXxb:" << dotXxb << endl;
    cout << "dotXxb/dotxb:" << dotXxb / dotxb << endl;
    end = MPI_Wtime();
    cout << "################" << endl;
    cout << "time cost: " << end - start << "s" << endl;
    MPI_Finalize();
    return 0;
}
/*
    LatticeFermi a(lat_x, lat_y, lat_z, lat_t, lat_s,lat_c);
    a.assign_random(686999);
    a.info();
    LatticeFermi b(lat_x, lat_y, lat_z, lat_t, lat_s,lat_c);
    b.info();
    b.assign_random(667);
    Complex dotXa;
    Complex dotXab;
    cout << "dotab:" << dot(a, b) << endl;
    cout << "dotXab:" << dotXab << endl;
    cout << "dotXab/dotab:" << dotXab / dot(a, b) << endl;
*/
/*
    LatticeGauge U(lat_x, lat_y, lat_z, lat_t, lat_s,lat_c);
    U.assign_random(0007);
    U.info();

    LatticePropagator prop = fermi_to_prop(a,2);
    print(0, 0, 0, 0, 1);
    print();
    prop.print(0, 0, 0, 1, 4);
    prop.print();
    U.print();
    LatticeFermi a1(lat_x, lat_y, lat_z, lat_t, lat_s,lat_c);
    a1.assign_random(12);
    LatticeFermi a2(lat_x, lat_y, lat_z, lat_t, lat_s,lat_c);
    a2.assign_random(13);
    LatticeFermi sum = a1 + a2;
    LatticeFermi diff = a1 - a2;
    LatticeFermi prod = a1 * a2;
    LatticeFermi quot = a1 / a2;
    LatticeFermi neg = -a1;
    a1.print(1, 1, 1, 1, 1);
    a2.print(1, 1, 1, 1, 1);
    // sum.print(1, 1, 1, 1, 1);
    // diff.print(1, 1, 1, 1, 1);
    // prod.print(1, 1, 1, 1, 1);
    // quot.print(1, 1, 1, 1, 1);
    // neg.print(1, 1, 1, 1, 1);

    bool same = (a1 == a2);
    bool different = (a1 != a2);
    cout << "same" << same;
    cout << "different" << different;

    a1.print();
    a2.print();
    cout << dot(a1, a2);
    cout << a1[0].imag;
    */
/*
LatticeFermi b(lat_x, lat_y, lat_z, lat_t, lat_s,lat_c);
// b.info();
b.assign_random(667);
int num_x(1), num_y(1), num_z(2), num_t(2), 1(1);
Complex dotXa;
Complex dotXab;
for (int x = 0; x < num_x; x++)
{
    for (int y = 0; y < num_y; y++)
    {
        for (int z = 0; z < num_z; z++)
        {
            for (int t = 0; t < num_t; t++)
            {
                for (int s = 0; s < 1; s++)
                {
                    int index = x * num_y * num_z * num_t * 1 + y * num_z * num_t * 1 + z * num_t * 1 + t * 1 + s;
                    if (rank == index)
                    {
                        LatticeFermi block_a = block(a, num_x, num_y, num_z, num_t, 1, x, y, z, t, s,c);
                        LatticeFermi block_b = block(b, num_x, num_y, num_z, num_t, 1, x, y, z, t, s,c);
                        dotXab =dotX(block_a, block_b);
                        cout << "block_dotab:" << dot(block_a,block_b) << endl;
                        dotXa =dotX(block_a, block_b);
                        cout << "block_dot:" << dot(block_a,block_b) << endl;

                    }
                }
            }
        }
    }
}
cout << "dotab:" << dot(a, b) << endl;
cout << "dotXab:" <<  dotXab << endl;
cout << "dotXab/dotab:" <<dotXab / dot(a, b) << endl;
// cout << "dot:" << dot(a, a) << endl;
// cout << "dotX:" <<  dotXa << endl;
// cout << "dotX/dot:" <<dotXa / dot(a, a) << endl;
*/
/*
LatticeFermi b0 = block(a, 1, 2, 1, 1, 1, 0, 0, 0, 0, 0);
LatticeFermi b1 = block(a, 1, 2, 1, 1, 1, 0, 1, 0, 0, 0);
cout << "dotx:" << dot(a, a) << endl;
cout << "norm2x:" << norm_2(a) << endl;
cout << "dotb0:" << dot(b0,b0) << endl;
cout << "dotb1:" << dot(b1,b1) << endl;
*/