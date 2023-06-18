#ifndef _LATTICE_VEC_HPP
#define _LATTICE_VEC_HPP
#include <complex.h>
#include <iostream>
// #include <cstdlib>

class lattice_vec
{
public:
    int size;
    lattice_vec(int size)
        : size(size)
    {
        _lattice_vec = new std::complex<double>[size];
    }
    ~lattice_vec()
    {
        delete [] _lattice_vec;
    }
    void clean()
    {
        fill(_lattice_vec, _lattice_vec + size, std::complex<double>(0, 0));
    }
    void clean1()
    {
        fill(_lattice_vec, _lattice_vec + size, std::complex<double>(1, 0));
    }
    void rand()
    {
        for (int i = 0; i < size; i++)
        {
            _lattice_vec[i] = std::complex<double>(std::rand() / 1e9, 0);
        }
    }
    std::complex<double> &operator[](int &i)
    {
        return _lattice_vec[i];
    }
    std::complex<double> &operator[](const int &i)
    {
        return _lattice_vec[i];
    }
    lattice_vec &operator=(lattice_vec &a)
    {
        copy(a._lattice_vec, a._lattice_vec + size, _lattice_vec);
        return *this;
    }
    lattice_vec &operator=(const lattice_vec &a)
    {
        copy(a._lattice_vec, a._lattice_vec + size, _lattice_vec);
        return *this;
    }
    lattice_vec operator+(lattice_vec &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] + a._lattice_vec[i];
        }
        return result;
    }
    lattice_vec operator+(const lattice_vec &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] + a._lattice_vec[i];
        }
        return result;
    }
    lattice_vec operator-(lattice_vec &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] - a._lattice_vec[i];
        }
        return result;
    }
    lattice_vec operator-(const lattice_vec &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] - a._lattice_vec[i];
        }
        return result;
    }
    lattice_vec operator*(double &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] * a;
        }
        return result;
    }
    lattice_vec operator*(const double &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] * a;
        }
        return result;
    }
    lattice_vec operator/(double &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] / a;
        }
        return result;
    }
    lattice_vec operator/(const double &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] / a;
        }
        return result;
    }
    lattice_vec operator*(std::complex<double> &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] * a;
        }
        return result;
    }
    lattice_vec operator*(const std::complex<double> &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] * a;
        }
        return result;
    }
    lattice_vec operator/(std::complex<double> &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] / a;
        }
        return result;
    }
    lattice_vec operator/(const std::complex<double> &a)
    {
        lattice_vec result(size);
        for (int i = 0; i < size; i++)
        {
            result[i] = this->_lattice_vec[i] / a;
        }
        return result;
    }
    std::complex<double> dot(lattice_vec &a)
    {
        std::complex<double> result(0.0, 0.0);
        for (int i = 0; i < size; i++)
        {
            result += a[i] * conj(this->_lattice_vec[i]);
        }
        return result;
    };
    double norm()
    {
        double result = 0.0;
        for (int i = 0; i < size; i++)
        {
            result += real(_lattice_vec[i] * conj(_lattice_vec[i]));
        }
        return result;
    }
    void print()
    {
        std::cout << "Lattice vector: " << _lattice_vec << std::endl;
        for (int i = 0; i < size; i++)
        {
            if (i > 24)
            {
                std::cout << "->...";
                break;
            }
            std::cout << _lattice_vec[i] << " ";
        }
        std::cout << std::endl;
    }

private:
    std::complex<double> *_lattice_vec;
};

#endif
