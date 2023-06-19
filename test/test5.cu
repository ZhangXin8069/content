#include <iostream>
#include <curand.h>
#include <curand_kernel.h>

// Define dimensions
const int DIM_X = 32;
const int DIM_Y = 32;
const int DIM_Z = 32;
const int DIM_T = 32;
const int NUM_PARITIES = 2;
const int VOLUME = DIM_X * DIM_Y * DIM_Z * DIM_T;
// Complex number structure
struct Complex
{
    double real;
    double imag;

    __device__ Complex operator*(const Complex &other) const
    {
        Complex result;
        result.real = real * other.real - imag * other.imag;
        result.imag = real * other.imag + imag * other.real;
        return result;
    }

    __device__ Complex operator+(const Complex &other) const
    {
        Complex result;
        result.real = real + other.real;
        result.imag = imag + other.imag;
        return result;
    }

    __device__ Complex operator-(const Complex &other) const
    {
        Complex result;
        result.real = real - other.real;
        result.imag = imag - other.imag;
        return result;
    }
};

__device__ Complex conj(const Complex &c)
{
    Complex result;
    result.real = c.real;
    result.imag = -c.imag;
    return result;
}

// Fermi field class
class FermiField
{
private:
    Complex *field;

public:
    int numParities;
    __host__ __device__ FermiField(Complex *fieldPtr, int parities) : field(fieldPtr), numParities(parities) {}

    __host__ __device__ Complex &getField(int index, int parity)
    {
        return field[parity * VOLUME + index];
    }

    __host__ __device__ void setField(int index, int parity, const Complex &value)
    {
        field[parity * VOLUME + index] = value;
    }
};

// Gauge field class
class GaugeField
{
private:
    Complex *field;
    int numParities;

public:
    __host__ __device__ GaugeField(Complex *fieldPtr, int parities) : field(fieldPtr), numParities(parities) {}

    __device__ Complex &getLink(int index, int mu, int parity)
    {
        return field[mu * (numParities * VOLUME) + parity * VOLUME + index];
    }

    __device__ void setLink(int index, int mu, int parity, const Complex &value)
    {
        field[mu * (numParities * VOLUME) + parity * VOLUME + index] = value;
    }
};

__global__ void setupRandomGenerator(curandState *devStates)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < VOLUME)
    {
        curand_init(1234, idx, 0, &devStates[idx]);
    }
}

__global__ void initInputFields(Complex *devFermiField, Complex *devGaugeField, curandState *devStates)
{
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < VOLUME)
    {
        int parity = idx % NUM_PARITIES;

        curandState localState = devStates[idx];

        devFermiField[parity * VOLUME + idx] = {curand_uniform(&localState), curand_uniform(&localState)};

        for (int dir = 0; dir < 4; dir++)
        {
            devGaugeField[(dir * NUM_PARITIES + parity) * VOLUME + idx] = {curand_uniform(&localState), curand_uniform(&localState)};
        }

        devStates[idx] = localState;
    }
}

// Kernel for the Dslash4 operation
__global__ void Dslash4(FermiField fermiField, GaugeField gaugeField, FermiField resultField, int parity)
{
    bool dag = true;
    const double a = 2.0;
    const double mass = 1.0;
    Complex i = {0.0, 1};
    Complex Half = {0.5, 0.0};
    Complex flag = {(dag == true) ? -1.0 : 1.0, 0};
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < VOLUME)
    {
        int x = idx % DIM_X;
        int y = (idx / DIM_X) % DIM_Y;
        int z = (idx / (DIM_X * DIM_Y)) % DIM_Z;
        int t = (idx / (DIM_X * DIM_Y * DIM_Z)) % DIM_T;

        // mass term
        for (int s = 0; s < fermiField.numParities; s++)
        {
            Complex value = fermiField.getField(idx, s);
            resultField.setField(idx, s, value *Complex{-(a + mass), 0.0});
        }

        // backward x
        int b_x = (x + DIM_X - 1) % DIM_X;
        Complex tmp = (fermiField.getField(idx, 0) + flag * fermiField.getField(idx, 1)) * Half * gaugeField.getLink(idx, 0, parity);
        resultField.setField(b_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 0, resultField.getField(b_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 0) + tmp);
        resultField.setField(b_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 1, resultField.getField(b_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 1) + flag * tmp);

        // forward x
        int f_x = (x + 1) % DIM_X;
        tmp = (fermiField.getField(idx, 0) - flag * fermiField.getField(idx, 1)) * Half * conj(gaugeField.getLink(f_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 0, parity));
        resultField.setField(f_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 0, resultField.getField(f_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 0) + tmp);
        resultField.setField(f_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 1, resultField.getField(f_x + DIM_X * (y + DIM_Y * (z + DIM_Z * t)), 1) - flag * tmp);

        // backward y
        int b_y = (y + DIM_Y - 1) % DIM_Y;
        tmp = (fermiField.getField(idx, 0) + flag * i * fermiField.getField(idx, 1)) * Half * gaugeField.getLink(idx, 1, parity);
        resultField.setField(x + DIM_X * (b_y + DIM_Y * (z + DIM_Z * t)), 0, resultField.getField(x + DIM_X * (b_y + DIM_Y * (z + DIM_Z * t)), 0) + tmp);
        resultField.setField(x + DIM_X * (b_y + DIM_Y * (z + DIM_Z * t)), 1, resultField.getField(x + DIM_X * (b_y + DIM_Y * (z + DIM_Z * t)), 1) - flag * i * tmp);

        // forward y
        int f_y = (y + 1) % DIM_Y;
        tmp = (fermiField.getField(idx, 0) - flag * i * fermiField.getField(idx, 1)) * Half * conj(gaugeField.getLink(x + DIM_X * (f_y + DIM_Y * (z + DIM_Z * t)), 1, parity));
        resultField.setField(x + DIM_X * (f_y + DIM_Y * (z + DIM_Z * t)), 0, resultField.getField(x + DIM_X * (f_y + DIM_Y * (z + DIM_Z * t)), 0) + tmp);
        resultField.setField(x + DIM_X * (f_y + DIM_Y * (z + DIM_Z * t)), 1, resultField.getField(x + DIM_X * (f_y + DIM_Y * (z + DIM_Z * t)), 1) + flag * i * tmp);

        // backward z
        int b_z = (z + DIM_Z - 1) % DIM_Z;
        tmp = (fermiField.getField(idx, 0) + flag * fermiField.getField(idx, 1)) * Half * gaugeField.getLink(idx, 2, parity);
        resultField.setField(x + DIM_X * (y + DIM_Y * (b_z + DIM_Z * t)), 0, resultField.getField(x + DIM_X * (y + DIM_Y * (b_z + DIM_Z * t)), 0) + tmp);
        resultField.setField(x + DIM_X * (y + DIM_Y * (b_z + DIM_Z * t)), 1, resultField.getField(x + DIM_X * (y + DIM_Y * (b_z + DIM_Z * t)), 1) + flag * tmp);

        // forward z
        int f_z = (z + 1) % DIM_Z;
        tmp = (fermiField.getField(idx, 0) - flag * fermiField.getField(idx, 1)) * Half * conj(gaugeField.getLink(x + DIM_X * (y + DIM_Y * (f_z + DIM_Z * t)), 2, parity));
        resultField.setField(x + DIM_X * (y + DIM_Y * (f_z + DIM_Z * t)), 0, resultField.getField(x + DIM_X * (y + DIM_Y * (f_z + DIM_Z * t)), 0) + tmp);
        resultField.setField(x + DIM_X * (y + DIM_Y * (f_z + DIM_Z * t)), 1, resultField.getField(x + DIM_X * (y + DIM_Y * (f_z + DIM_Z * t)), 1) - flag * tmp);

        // backward t
        int b_t = (t + DIM_T - 1) % DIM_T;
        tmp = (fermiField.getField(idx, 0) + flag * i * fermiField.getField(idx, 1)) * Half * gaugeField.getLink(idx, 3, parity);
        resultField.setField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * b_t)), 0, resultField.getField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * b_t)), 0) + tmp);
        resultField.setField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * b_t)), 1, resultField.getField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * b_t)), 1) - flag * i * tmp);

        // forward t
        int f_t = (t + 1) % DIM_T;
        tmp = (fermiField.getField(idx, 0) - flag * i * fermiField.getField(idx, 1)) * Half * conj(gaugeField.getLink(x + DIM_X * (y + DIM_Y * (z + DIM_Z * f_t)), 3, parity));
        resultField.setField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * f_t)), 0, resultField.getField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * f_t)), 0) + tmp);
        resultField.setField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * f_t)), 1, resultField.getField(x + DIM_X * (y + DIM_Y * (z + DIM_Z * f_t)), 1) + flag * i * tmp);
    }
}

int main()
{
    // Allocate memory on the host
    Complex *hostFermiField = new Complex[NUM_PARITIES * VOLUME];
    Complex *hostGaugeField = new Complex[4 * NUM_PARITIES * VOLUME];
    curandState *hostRandomStates = new curandState[VOLUME];

    // Allocate memory on the device
    Complex *devFermiField;
    Complex *devGaugeField;
    curandState *devRandomStates;
    cudaMalloc(&devFermiField, sizeof(Complex) * NUM_PARITIES * VOLUME);
    cudaMalloc(&devGaugeField, sizeof(Complex) * 4 * NUM_PARITIES * VOLUME);
    cudaMalloc(&devRandomStates, sizeof(curandState) * VOLUME);

    // Initialize random number generator states on the device
    int numThreads = 256;
    int numBlocks = (VOLUME + numThreads - 1) / numThreads;
    setupRandomGenerator<<<numBlocks, numThreads>>>(devRandomStates);
    cudaDeviceSynchronize();

    // Initialize input fields on the device
    initInputFields<<<numBlocks, numThreads>>>(devFermiField, devGaugeField, devRandomStates);
    cudaDeviceSynchronize();

    // Copy input fields from device to host
    cudaMemcpy(hostFermiField, devFermiField, sizeof(Complex) * NUM_PARITIES * VOLUME, cudaMemcpyDeviceToHost);
    cudaMemcpy(hostGaugeField, devGaugeField, sizeof(Complex) * 4 * NUM_PARITIES * VOLUME, cudaMemcpyDeviceToHost);

    // Perform Dslash4 operation on the device
    FermiField fermiField(devFermiField, NUM_PARITIES);
    GaugeField gaugeField(devGaugeField, NUM_PARITIES);
    FermiField resultField(devFermiField, NUM_PARITIES);

    Dslash4<<<numBlocks, numThreads>>>(fermiField, gaugeField, resultField, 0);
    cudaDeviceSynchronize();

    // Copy result field from device to host
    cudaMemcpy(hostFermiField, devFermiField, sizeof(Complex) * NUM_PARITIES * VOLUME, cudaMemcpyDeviceToHost);

    // Output results
    for (int p = 0; p < NUM_PARITIES; p++)
    {
        for (int x = 0; x < DIM_X; x++)
        {
            for (int y = 0; y < DIM_Y; y++)
            {
                for (int z = 0; z < DIM_Z; z++)
                {
                    for (int t = 0; t < DIM_T; t++)
                    {
                        Complex value = hostFermiField[p * VOLUME + x + DIM_X * (y + DIM_Y * (z + DIM_Z * t))];
                        std::cout << "Result[" << p << "][" << x << "][" << y << "][" << z << "][" << t
                                  << "]: (" << value.real << ", " << value.imag << ")" << std::endl;
                    }
                }
            }
        }
    }

    // Free memory on the device
    cudaFree(devFermiField);
    cudaFree(devGaugeField);
    cudaFree(devRandomStates);

    // Free memory on the host
    delete[] hostFermiField;
    delete[] hostGaugeField;
    delete[] hostRandomStates;

    return 0;
}
