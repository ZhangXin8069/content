#ifndef _ZVEXTOR_H
#define _ZVEXTOR_H
#include <iostream>
#include <complex>
#include <random>
#include <mpi.h>

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
#endif
