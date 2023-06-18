#ifndef _MY_VEXTOR_H
#define _MY_VEXTOR_H
#include <vector>
#include <complex>
#include <mpi.h>

class ComplexVector
{
public:
    explicit ComplexVector(int size) : data_(size)
    {
        this->clean_0();
        MPI_Comm_rank(MPI_COMM_WORLD, &node_rank_);
        MPI_Comm_size(MPI_COMM_WORLD, &node_size_);
    }

    ComplexVector(const ComplexVector &other) : data_(other.data_) {}

    // Give size
    int size() const { return data_.size(); }

    // Turn into a other ComplexVector
    void turn_into(const ComplexVector &other)
    {
        data_.resize(other.size());
        for (int i = 0; i < data_.size(); i++)
        {
            data_[i] = other[i];
        }
    }

    // Clean into 0
    void clean_0()
    {
        for (int i = 0; i < data_.size(); i++)
        {
            data_[i] = std::complex<double>(0, 0);
        }
    }

    // Clean into 1
    void clean_1()
    {
        for (int i = 0; i < data_.size(); i++)
        {
            data_[i] = std::complex<double>(1, 0);
        }
    }

    // Clean into rand
    void clean_rand()
    {
        for (int i = 0; i < data_.size(); i++)
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
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] + other[i];
        }
        return result;
    }

    // Add ComplexVector with a complex number
    ComplexVector operator+(const std::complex<double> &scalar) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] + scalar;
        }
        return result;
    }

    // Subtract two ComplexVectors
    ComplexVector operator-(const ComplexVector &other) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] - other[i];
        }
        return result;
    }

    // Subtract ComplexVector with a complex number
    ComplexVector operator-(const std::complex<double> &scalar) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] - scalar;
        }
        return result;
    }

    // Multiply two ComplexVectors element-wise
    ComplexVector operator*(const ComplexVector &other) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] * other[i];
        }
        return result;
    }

    // Multiply ComplexVector with a complex number
    ComplexVector operator*(const std::complex<double> &scalar) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] * scalar;
        }
        return result;
    }

    // Divide two ComplexVectors element-wise
    ComplexVector operator/(const ComplexVector &other) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] / other[i];
        }
        return result;
    }

    // Divide ComplexVector with a complex number
    ComplexVector operator/(const std::complex<double> &scalar) const
    {
        ComplexVector result(data_.size());
        for (int i = 0; i < data_.size(); i++)
        {
            result[i] = data_[i] / scalar;
        }
        return result;
    }

    // Calculate dot product not using OpenMPI non-blocking communication
    std::complex<double> dot(const ComplexVector &other)
    {
        std::complex<double> result(0.0, 0.0);
        for (int i = 0; i < data_.size(); i++)
        {
            result += data_[i] * conj(other[i]);
        }
        return result;
    };

    // Calculate dot product using OpenMPI non-blocking communication
    std::complex<double> dotX(const ComplexVector &other) const
    {
        std::complex<double> tmp(0.0, 0.0);
        for (int i = 0; i < data_.size(); i++)
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
        for (int i = 0; i < data_.size(); i++)
        {
            result += data_[i] * conj(data_[i]);
        }
        return result;
    };

    // Calculate norm2 product using OpenMPI non-blocking communication
    std::complex<double> norm2X() const
    {
        std::complex<double> tmp(0.0, 0.0);
        for (int i = 0; i < data_.size(); i++)
        {
            tmp += data_[i] * conj(data_[i]);
        }
        MPI_Barrier(MPI_COMM_WORLD);
        MPI_Allreduce(&tmp, &tmp, 1, MPI_DOUBLE_COMPLEX, MPI_SUM, MPI_COMM_WORLD);

        return tmp;
    }

    // private:
    std::vector<std::complex<double>> data_;
    int node_rank_, node_size_;
};

class Vector
{
public:
    Vector(int size) : size(size), data(size) {}
    Vector(const std::vector<double> &v) : size(v.size()), data(v) {}

    int size;
    std::vector<double> data;

    double operator[](int i) const { return data[i]; }
    double &operator[](int i) { return data[i]; }

    Vector operator+(const Vector &other) const
    {
        Vector result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = data[i] + other[i];
        }
        return result;
    }

    Vector operator-(const Vector &other) const
    {
        Vector result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = data[i] - other[i];
        }
        return result;
    }

    Vector operator*(double scalar) const
    {
        Vector result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = data[i] * scalar;
        }
        return result;
    }

    double dot(const Vector &other) const
    {
        double result = 0;
        for (int i = 0; i < size; i++)
        {
            result += data[i] * other[i];
        }
        return result;
    }
};

std::ostream &operator<<(std::ostream &stream, const Vector &v)
{
    stream << "[";
    for (int i = 0; i < v.size; i++)
    {
        stream << v[i];
        if (i > 4)
        {
            stream << "->...";
            break;
        }
        if (i != v.size - 1)
        {
            stream << ", ";
        }
    }
    stream << "]\n";
    return stream;
}
std::ostream &operator<<(std::ostream &stream, const ComplexVector &v)
{
    stream << "[";
    for (int i = 0; i < v.size(); i++)
    {
        stream << v[i];
        if (i > 4)
        {
            stream << "->...";
            break;
        }
        if (i != v.size() - 1)
        {
            stream << ", ";
        }
    }
    stream << "]\n";
    return stream;
}
#endif
