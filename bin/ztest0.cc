#include <iostream>
#include <random>
#include <mpi.h>
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
    LatticeFermi(const int &lat_x, const int &lat_y, const int &lat_z, const int &lat_t, const int &lat_s, const int &lat_c)
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
    const Complex &operator[](const int &index) const
    {
        return lattice_vec[index];
    }
    Complex &operator[](const int &index)
    {
        return lattice_vec[index];
    }
    const Complex &operator()(const int &index_x, const int &index_y, const int &index_z, const int &index_t, const int &index_s, const int &index_c) const
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c + index_y * lat_z * lat_t * lat_s * lat_c + index_z * lat_t * lat_s * lat_c + index_t * lat_s * lat_c + index_s * lat_c + index_c;
        return lattice_vec[index];
    }
    Complex &operator()(const int &index_x, const int &index_y, const int &index_z, const int &index_t, const int &index_s, const int &index_c)
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
    void print(const int &index_x, const int &index_y, const int &index_z, const int &index_t, const int &index_s, const int &index_c)
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
    LatticeFermi block(const int &num_x, const int &num_y, const int &num_z, const int &num_t, const int &index_x, const int &index_y, const int &index_z, const int &index_t)
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
    LatticeFermi block(const int &num_x, const int &num_y, const int &num_z, const int &num_t)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return block(num_x, num_y, num_z, num_t,
                     rank / (num_y * num_z * num_t),
                     (rank / (num_z * num_t)) % num_y,
                     (rank / num_t) % num_z,
                     rank % num_t);
    }
    LatticeFermi reback(const int &num_x, const int &num_y, const int &num_z, const int &num_t)
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
    LatticeGauge(const int &lat_x, const int &lat_y, const int &lat_z, const int &lat_t, const int &lat_s, const int &lat_c)
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
    const Complex &operator[](const int &index) const
    {
        return lattice_vec[index];
    }
    Complex &operator[](const int &index)
    {
        return lattice_vec[index];
    }
    const Complex &operator()(const int &index_x, const int &index_y, const int &index_z, const int &index_t, const int &index_s, const int &index_c0, const int &index_c1) const
    {
        int index = index_x * lat_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_y * lat_z * lat_t * lat_s * lat_c0 * lat_c1 + index_z * lat_t * lat_s * lat_c0 * lat_c1 + index_t * lat_s * lat_c0 * lat_c1 + index_s * lat_c0 * lat_c1 + index_c0 * lat_c1 + index_c1;
        return lattice_vec[index];
    }
    Complex &operator()(const int &index_x, const int &index_y, const int &index_z, const int &index_t, const int &index_s, const int &index_c0, const int &index_c1)
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
    void print(const int &index_x, const int &index_y, const int &index_z, const int &index_t, const int &index_s, const int &index_c0, const int &index_c1)
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
    Complex dot(const LatticeGauge &other)
    {
        Complex result;
        for (int i = 0; i < size; i++)
        {
            result = result + lattice_vec[i].conj() * other[i];
        }
        return result;
    }
    Complex dotX(const LatticeGauge &other)
    {
        Complex local_result;
        Complex global_result;
        local_result = dot(other);
        MPI_Allreduce(&local_result, &global_result, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

        return global_result;
    }
    LatticeGauge block(const int &num_x, const int &num_y, const int &num_z, const int &num_t, const int &index_x, const int &index_y, const int &index_z, const int &index_t)
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
    LatticeGauge block(const int &num_x, const int &num_y, const int &num_z, const int &num_t)
    {
        int rank;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        return block(num_x, num_y, num_z, num_t,
                     rank / (num_y * num_z * num_t),
                     (rank / (num_z * num_t)) % num_y,
                     (rank / num_t) % num_z,
                     rank % num_t);
    }
    LatticeGauge reback(const int &num_x, const int &num_y, const int &num_z, const int &num_t)
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
void dslash(LatticeGauge &U, LatticeFermi &src, LatticeFermi &dest, const int &num_x, const int &num_y, const int &num_z, const int &num_t, const double mass = 1.0, const bool dag = false)
{
    dest.assign_zero();
    const double a = 2.0;
    const Complex i(0.0, 1.0);
    Complex tmp0[3];
    Complex tmp1[3];
    Complex g0[2];
    Complex g1[2];
    int s0[2];
    int s1[2];
    int d;
    const double Half = 0.5;
    double flag = (dag == true) ? -1 : 1;
    Complex flag0;
    Complex flag1;
    for (int x = 0; x < U.lat_x; x++)
    {
        for (int y = 0; y < U.lat_y; y++)
        {
            for (int z = 0; z < U.lat_z; z++)
            {
                for (int t = 0; t < U.lat_t; t++)
                {
                    // // mass term and others
                    // for (int s = 0; s < U.lat_s; s++)
                    // {
                    //     dest(x, y, z, t, s, c) += -src(x, y, z, t, s, c) * (a + mass);
                    // }

                    // backward x
                    int b_x = (x + U.lat_x - 1) % U.lat_x;
                    d = 0;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 3;
                    g0[1] = i;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 2;
                    g1[1] = i;
                    flag0 = -i;
                    flag1 = -i;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(b_x, y, z, t, s0[0], c1) * g0[0] + src(b_x, y, z, t, s0[1], c1) * g0[1]) * U(b_x, y, z, t, d, c0, c1).conj();
                            tmp1[c0] += (src(b_x, y, z, t, s1[0], c1) * g1[0] + src(b_x, y, z, t, s1[1], c1) * g1[1]) * U(b_x, y, z, t, d, c0, c1).conj();
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp1[c0] * flag1;
                        dest(x, y, z, t, 3, c0) += tmp0[c0] * flag0;
                    }
                    // forward x
                    int f_x = (x + 1) % U.lat_x;
                    d = 0;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 3;
                    g0[1] = -i;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 2;
                    g1[1] = -i;
                    flag0 = i;
                    flag1 = i;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(f_x, y, z, t, s0[0], c1) * g0[0] + src(f_x, y, z, t, s0[1], c1) * g0[1]) * U(x, y, z, t, d, c0, c1);
                            tmp1[c0] += (src(f_x, y, z, t, s1[0], c1) * g1[0] + src(f_x, y, z, t, s1[1], c1) * g1[1]) * U(x, y, z, t, d, c0, c1);
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp1[c0] * flag1;
                        dest(x, y, z, t, 3, c0) += tmp0[c0] * flag0;
                    }
                    // backward y
                    int b_y = (y + U.lat_y - 1) % U.lat_y;
                    d = 1;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 3;
                    g0[1] = -1;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 2;
                    g1[1] = 1;
                    flag0 = -1;
                    flag1 = 1;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(x, b_y, z, t, s0[0], c1) * g0[0] + src(x, b_y, z, t, s0[1], c1) * g0[1]) * U(x, b_y, z, t, d, c0, c1).conj();
                            tmp1[c0] += (src(x, b_y, z, t, s1[0], c1) * g1[0] + src(x, b_y, z, t, s1[1], c1) * g1[1]) * U(x, b_y, z, t, d, c0, c1).conj();
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp1[c0] * flag1;
                        dest(x, y, z, t, 3, c0) += tmp0[c0] * flag0;
                    }
                    // forward y
                    int f_y = (y + 1) % U.lat_y;
                    d = 1;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 3;
                    g0[1] = 1;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 2;
                    g1[1] = -1;
                    flag0 = 1;
                    flag1 = -1;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(x, f_y, z, t, s0[0], c1) * g0[0] + src(x, f_y, z, t, s0[1], c1) * g0[1]) * U(x, y, z, t, d, c0, c1);
                            tmp1[c0] += (src(x, f_y, z, t, s1[0], c1) * g1[0] + src(x, f_y, z, t, s1[1], c1) * g1[1]) * U(x, y, z, t, d, c0, c1);
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp1[c0] * flag1;
                        dest(x, y, z, t, 3, c0) += tmp0[c0] * flag0;
                    }
                    // backward z
                    int b_z = (z + U.lat_z - 1) % U.lat_z;
                    d = 2;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 2;
                    g0[1] = i;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 3;
                    g1[1] = -i;
                    flag0 = -i;
                    flag1 = i;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(x, y, b_z, t, s0[0], c1) * g0[0] + src(x, y, b_z, t, s0[1], c1) * g0[1]) * U(x, y, b_z, t, d, c0, c1).conj();
                            tmp1[c0] += (src(x, y, b_z, t, s1[0], c1) * g1[0] + src(x, y, b_z, t, s1[1], c1) * g1[1]) * U(x, y, b_z, t, d, c0, c1).conj();
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp0[c0] * flag0;
                        dest(x, y, z, t, 3, c0) += tmp1[c0] * flag1;
                    }
                    // forward z
                    int f_z = (z + 1) % U.lat_z;
                    d = 2;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 2;
                    g0[1] = -i;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 3;
                    g1[1] = i;
                    flag0 = i;
                    flag1 = -i;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(x, y, f_z, t, s0[0], c1) * g0[0] + src(x, y, f_z, t, s0[1], c1) * g0[1]) * U(x, y, z, t, d, c0, c1);
                            tmp1[c0] += (src(x, y, f_z, t, s1[0], c1) * g1[0] + src(x, y, f_z, t, s1[1], c1) * g1[1]) * U(x, y, z, t, d, c0, c1);
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp0[c0] * flag0;
                        dest(x, y, z, t, 3, c0) += tmp1[c0] * flag1;
                    }
                    // backward t
                    int b_t = (t + U.lat_t - 1) % U.lat_t;
                    d = 3;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 2;
                    g0[1] = 1;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 3;
                    g1[1] = 1;
                    flag0 = 1;
                    flag1 = 1;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(x, y, z, b_t, s0[0], c1) * g0[0] + src(x, y, z, b_t, s0[1], c1) * g0[1]) * U(x, y, z, b_t, d, c0, c1).conj();
                            tmp1[c0] += (src(x, y, z, b_t, s1[0], c1) * g1[0] + src(x, y, z, b_t, s1[1], c1) * g1[1]) * U(x, y, z, b_t, d, c0, c1).conj();
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp0[c0] * flag0;
                        dest(x, y, z, t, 3, c0) += tmp1[c0] * flag1;
                    }
                    // forward t
                    int f_t = (t + 1) % U.lat_t;
                    d = 3;
                    tmp0[0] = 0;
                    tmp0[1] = 0;
                    tmp0[2] = 0;
                    tmp1[0] = 0;
                    tmp1[1] = 0;
                    tmp1[2] = 0;
                    s0[0] = 0;
                    g0[0] = 1;
                    s0[1] = 2;
                    g0[1] = -1;
                    s1[0] = 1;
                    g1[0] = 1;
                    s1[1] = 3;
                    g1[1] = -1;
                    flag0 = -1;
                    flag1 = -1;
                    for (int c0 = 0; c0 < U.lat_c0; c0++)
                    {
                        for (int c1 = 0; c1 < U.lat_c1; c1++)
                        {
                            tmp0[c0] += (src(x, y, z, f_t, s0[0], c1) * g0[0] + src(x, y, z, f_t, s0[1], c1) * g0[1]) * U(x, y, z, t, d, c0, c1);
                            tmp1[c0] += (src(x, y, z, f_t, s1[0], c1) * g1[0] + src(x, y, z, f_t, s1[1], c1) * g1[1]) * U(x, y, z, t, d, c0, c1);
                        }
                        dest(x, y, z, t, 0, c0) += tmp0[c0];
                        dest(x, y, z, t, 1, c0) += tmp1[c0];
                        dest(x, y, z, t, 2, c0) += tmp0[c0] * flag0;
                        dest(x, y, z, t, 3, c0) += tmp1[c0] * flag1;
                    }
                }
            }
        }
    }
}
void dslash(LatticeGauge &U, LatticeFermi &src, LatticeFermi &dest, const double mass = 1.0, const bool dag = true)
{
    dest = src * 500 + 0.3;
}
void cg(LatticeGauge &U0, LatticeFermi &b0, LatticeFermi &x0, const int &num_x, const int &num_y, const int &num_z, const int &num_t, const double mass = 1.0, const bool dag = true, const int MAX_ITER = 1e6, double TOL = 1e-6)
{
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Complex rho_prev(1.0, 0.0), rho(0.0, 0.0), alpha(1.0, 0.0), omega(1.0, 0.0), beta(0.0, 0.0);
    double r_norm2 = 0;
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
        dslash(U, p, v, mass, dag);
        alpha = rho / r_tilde.dotX(v);
        s = r - v * alpha;
        // t = A * s;
        // dslash(U, s, t, num_x, num_y, num_z, num_t, mass, dag);
        dslash(U, s, t, mass, dag);
        omega = t.dotX(s) / t.dotX(t);
        x = x + p * alpha + s * omega;
        r = s - t * omega;
        r_norm2 = r.norm_2X();
        std::cout << "##loop "
                  << i
                  << "##Residual:"
                  << r_norm2
                  << std::endl;
        if (r_norm2 < TOL || i == MAX_ITER - 1)
        {
            x0 = x.reback(num_x, num_y, num_z, num_t);
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
    int lat_x(8), lat_y(8), lat_z(8), lat_t(32), lat_s(4), lat_c(3);
    int num_x(1), num_y(1), num_z(2), num_t(4);
    int MAX_ITER(1e6);
    double TOL(1e-6);
    LatticeGauge U(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    LatticeFermi b(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    LatticeFermi x(lat_x, lat_y, lat_z, lat_t, lat_s, lat_c);
    U.assign_random(111);
    b.assign_random(222);
    x.assign_random(333);
    cg(U, b, x, num_x, num_y, num_z, num_t);
    x = x.reback(num_x, num_y, num_z, num_t);
    std::cout << "######x.norm_2():" << x.norm_2() << " ######" << std::endl;
    // if(rank==0){
    // x.print();
    // }
    end = MPI_Wtime();
    std::cout << "######time cost:" << end - start << "s ######" << std::endl;
    MPI_Finalize();
    return 0;
}
