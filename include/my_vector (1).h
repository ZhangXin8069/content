#ifndef _MY_VEXTOR_H
#define _MY_VEXTOR_H
#include <vector>
#include <complex>
#include <mpi.h>

class ComplexVector
{
public:
    explicit ComplexVector(int size) : size_(size)
    {
        data_ = new std::complex<double>[size_];
        MPI_Comm_rank(MPI_COMM_WORLD, &node_rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &node_size_);
    }

    ComplexVector(const ComplexVector &other) : data_(other.data_)
    {
    }

    // Clean into 0
    void clean_0()
    {
        for (int i = 0; i < size_; i++)
        {
            data_[i] = std::complex<double>(0, 0);
        }
    }

    // Clean into 1
    void clean_1()
    {
        for (int i = 0; i < size_; i++)
        {
            data_[i] = std::complex<double>(1, 0);
        }
    }

    // Clean into rand
    void clean_rand()
    {
        for (int i = 0; i < size_; i++)
        {
            data_[i] = std::complex<double>(rand() / double(RAND_MAX), rand() / double(RAND_MAX));
        }
    }

    // Access elements of the vector using [] operator
    std::complex<double> operator[](int index) const
    {
        return data_[index];
    }
    std::complex<double> &operator[](int index)
    {
        return data_[index];
    }

    // Add two ComplexVectors
    ComplexVector operator+(const ComplexVector &other) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] + other[i];
        }
        return result;
    }

    // Add ComplexVector with a complex number
    ComplexVector operator+(const std::complex<double> &scalar) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] + scalar;
        }
        return result;
    }

    // Subtract two ComplexVectors
    ComplexVector operator-(const ComplexVector &other) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] - other[i];
        }
        return result;
    }

    // Subtract ComplexVector with a complex number
    ComplexVector operator-(const std::complex<double> &scalar) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] - scalar;
        }
        return result;
    }

    // Multiply two ComplexVectors element-wise
    ComplexVector operator*(const ComplexVector &other) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] * other[i];
        }
        return result;
    }

    // Multiply ComplexVector with a complex number
    ComplexVector operator*(const std::complex<double> &scalar) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] * scalar;
        }
        return result;
    }

    // Divide two ComplexVectors element-wise
    ComplexVector operator/(const ComplexVector &other) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] / other[i];
        }
        return result;
    }

    // Divide ComplexVector with a complex number
    ComplexVector operator/(const std::complex<double> &scalar) const
    {
        ComplexVector result(size_);
        for (int i = 0; i < size_; i++)
        {
            result[i] = data_[i] / scalar;
        }
        return result;
    }

    // Calculate dot product not using OpenMPI non-blocking communication
    std::complex<double> dot(const ComplexVector &other)
    {
        std::complex<double> result(0.0, 0.0);
        for (int i = 0; i < size_; i++)
        {
            result += data_[i] * conj(other[i]);
        }
        return result;
    };

    // Calculate dot product using OpenMPI non-blocking communication
    std::complex<double> dotX(const ComplexVector &other) const
    {
        std::complex<double> tmp(0.0, 0.0);
        for (int i = 0; i < size_; i++)
        {
            tmp += data_[i] * conj(other[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmp, &tmp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

        return tmp;
    }

    // Calculate norm2 product not using OpenMPI non-blocking communication
    std::complex<double> norm2()
    {
        std::complex<double> result(0.0, 0.0);
        for (int i = 0; i < size_; i++)
        {
            result += data_[i] * conj(data_[i]);
        }
        return result;
    };

    // Calculate norm2 product using OpenMPI non-blocking communication
    std::complex<double> norm2X() const
    {
        std::complex<double> tmp(0.0, 0.0);
        for (int i = 0; i < size_; i++)
        {
            tmp += data_[i] * conj(data_[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmp, &tmp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

        return tmp;
    }

    // private:
    std::complex<double> *data_;
    int size_, node_rank_, node_size_;
};

std::ostream &operator<<(std::ostream &stream, const ComplexVector &v)
{
    stream << "[";
    for (int i = 0; i < v.size_; i++)
    {
        stream << v[i];
        if (i > 4)
        {
            stream << "->...";
            break;
        }
        if (i != v.size_ - 1)
        {
            stream << ", ";
        }
    }
    stream << "]\n";
    return stream;
}
#endif
