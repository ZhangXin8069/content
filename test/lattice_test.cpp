#include "lattice_vec.hpp"
int main()
{
    const int size = 25;
    lattice_vec vec1(size);
    lattice_vec vec2(size);
    for (int i = 0; i < vec1.size; i++)
        vec1[i] = std::complex<double>(i, i + 1);
    for (int i = 0; i < vec2.size; i++)
        vec2[i] = std::complex<double>(i, i + 2);
    std::cout << "Norm: " << vec1.norm() << std::endl;
    std::cout << "Norm: " << vec2.norm() << std::endl;
    vec1.print();
    vec2.print();
    // vec1 = vec2;
    vec1.print();
    vec2.print();
    lattice_vec vec3 = vec1 + vec2;
    lattice_vec vec4 = vec1 - vec2;
    vec3.print();
    vec4.print();
    lattice_vec vec5 = vec1 / 2;
    lattice_vec vec6 = vec1 * 2;
    vec5.print();
    vec6.print();

    std::cout << "Dot product: " << vec1.dot(vec2) << std::endl;
    vec1.clean();
    vec1.rand();
    vec1.print();
    return 0;
}
