#include <iostream>
#include <vector>
#include <complex>

class ComplexVector
{
public:
    explicit ComplexVector(int size) : data_(size) {}

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

    int size() const { return data_.size(); }

private:
    std::vector<std::complex<double>> data_;
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
        if (i != v.size - 1)
        {
            stream << ", ";
        }
    }
    stream << "]";
    return stream;
}
