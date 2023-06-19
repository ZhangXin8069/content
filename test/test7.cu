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
    float real;
    float imag;

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

    // You can also define other arithmetic operators if needed
};

// Fermi field class
class FermiField
{
private:
    Complex *field;
    int numParities;

public:
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

// Kernel for the dslash operation
__global__ void dslash(FermiField fermiField, GaugeField gaugeField, FermiField resultField, int parity)
{
    int index = blockIdx.x * blockDim.x + threadIdx.x;

    if (index < VOLUME)
    {
        int x = index % DIM_X;
        int y = (index / DIM_X) % DIM_Y;
        int z = (index / (DIM_X * DIM_Y)) % DIM_Z;
        int t = (index / (DIM_X * DIM_Y * DIM_Z)) % DIM_T;

        Complex result = {0.0f, 0.0f};

        for (int mu = 0; mu < 4; mu++)
        {
            int forwardIndex = index;
            int backwardIndex = index;

            // Forward direction
            if (x == DIM_X - 1 && mu == 0)
                forwardIndex = index + 1;
            else if (x < DIM_X - 1 && mu == 0)
                forwardIndex = index + DIM_X;
            else if (y == DIM_Y - 1 && mu == 1)
                forwardIndex = index + DIM_X * DIM_Y;
            else if (y < DIM_Y - 1 && mu == 1)
                forwardIndex = index + DIM_X;
            else if (z == DIM_Z - 1 && mu == 2)
                forwardIndex = index + DIM_X * DIM_Y * DIM_Z;
            else if (z < DIM_Z - 1 && mu == 2)
                forwardIndex = index + DIM_X * DIM_Y * DIM_Z;
            else if (t == DIM_T - 1 && mu == 3)
                forwardIndex = index + DIM_X * DIM_Y * DIM_Z * DIM_T;
            else if (t < DIM_T - 1 && mu == 3)
                forwardIndex = index + DIM_X * DIM_Y * DIM_Z * DIM_T;

            result = result + gaugeField.getLink(index, mu, parity) * fermiField.getField(forwardIndex, 1 - parity);

            // Backward direction
            if (x == 0 && mu == 0)
                backwardIndex = index - 1;
            else if (x > 0 && mu == 0)
                backwardIndex = index - DIM_X;
            else if (y == 0 && mu == 1)
                backwardIndex = index - DIM_X * DIM_Y;
            else if (y > 0 && mu == 1)
                backwardIndex = index - DIM_X;
            else if (z == 0 && mu == 2)
                backwardIndex = index - DIM_X * DIM_Y * DIM_Z;
            else if (z > 0 && mu == 2)
                backwardIndex = index - DIM_X * DIM_Y * DIM_Z;
            else if (t == 0 && mu == 3)
                backwardIndex = index - DIM_X * DIM_Y * DIM_Z * DIM_T;
            else if (t > 0 && mu == 3)
                backwardIndex = index - DIM_X * DIM_Y * DIM_Z * DIM_T;

            result = result + gaugeField.getLink(backwardIndex, mu, 1 - parity) * fermiField.getField(backwardIndex, 1 - parity);
        }

        resultField.setField(index, parity, result);
    }
}

// Stable Conjugate Gradient algorithm
void conjugateGradient(FermiField fermiField, GaugeField gaugeField, FermiField resultField, int maxIterations, float tolerance)
{
    FermiField residualField = resultField;

    for (int iteration = 0; iteration < maxIterations; iteration++)
    {
        // Compute residual vector: r = b - A * x
        dslash<<<(VOLUME + 255) / 256, 256>>>(fermiField, gaugeField, resultField, 0);
        dslash<<<(VOLUME + 255) / 256, 256>>>(fermiField, gaugeField, resultField, 1);

        // Update residual vector: r = b - A * x
        // (Note: Here, we assume x = 0, so r = b)

        // Compute squared norm of the residual: r_norm = |r|^2
        float rNorm = 0.0f;
        for (int index = 0; index < VOLUME; index++)
        {
            Complex residual = residualField.getField(index, 0);
            rNorm += residual.real * residual.real + residual.imag * residual.imag;
        }

        // Check convergence
        if (rNorm < tolerance)
        {
            break;
        }

        // Compute Ap
        dslash<<<(VOLUME + 255) / 256, 256>>>(fermiField, gaugeField, resultField, 0);
        dslash<<<(VOLUME + 255) / 256, 256>>>(fermiField, gaugeField, resultField, 1);

        // Compute squared norm of Ap: Ap_norm = |Ap|^2
        float apNorm = 0.0f;
        for (int index = 0; index < VOLUME; index++)
        {
            Complex ap = resultField.getField(index, 0);
            apNorm += ap.real * ap.real + ap.imag * ap.imag;
        }

        // Compute alpha: alpha = r_norm / Ap_norm
        float alpha = rNorm / apNorm;

        // Update solution vector: x = x + alpha * p
        for (int index = 0; index < VOLUME; index++)
        {
            Complex result = resultField.getField(index, 0);
            Complex residual = residualField.getField(index, 0);
            Complex solution = result;
            solution.real += alpha * residual.real;
            solution.imag += alpha * residual.imag;
            resultField.setField(index, 0, solution);
        }

        // Compute new residual vector: r = r - alpha * Ap
        for (int index = 0; index < VOLUME; index++)
        {
            Complex ap = resultField.getField(index, 0);
            Complex residual = residualField.getField(index, 0);
            Complex updatedResidual = residual;
            updatedResidual.real -= alpha * ap.real;
            updatedResidual.imag -= alpha * ap.imag;
            residualField.setField(index, 0, updatedResidual);
        }

        // Compute squared norm of the updated residual: r_norm_new = |r|^2
        float rNormNew = 0.0f;
        for (int index = 0; index < VOLUME; index++)
        {
            Complex updatedResidual = residualField.getField(index, 0);
            rNormNew += updatedResidual.real * updatedResidual.real + updatedResidual.imag * updatedResidual.imag;
        }

        // Check convergence
        if (rNormNew < tolerance)
        {
            break;
        }

        // Compute beta: beta = r_norm_new / r_norm
        float beta = rNormNew / rNorm;

        // Update direction vector: p = r + beta * p
        for (int index = 0; index < VOLUME; index++)
        {
            Complex residual = residualField.getField(index, 0);
            Complex direction = resultField.getField(index, 0);
            direction.real = residual.real + beta * direction.real;
            direction.imag = residual.imag + beta * direction.imag;
            resultField.setField(index, 0, direction);
        }

        // Update residual norm
        rNorm = rNormNew;
    }
}

int main()
{
    cudaDeviceSynchronize();
    // Random number generation setup
    curandState *devStates;
    cudaMalloc((void **)&devStates, VOLUME * sizeof(curandState));
    setupRandomGenerator<<<(VOLUME + 255) / 256, 256>>>(devStates);

    // Initialize input fields
    Complex *devFermiField;
    Complex *devGaugeField;
    Complex *devResultField;
    cudaMalloc((void **)&devFermiField, NUM_PARITIES * VOLUME * sizeof(Complex));
    cudaMalloc((void **)&devGaugeField, 4 * NUM_PARITIES * VOLUME * sizeof(Complex));
    cudaMalloc((void **)&devResultField, NUM_PARITIES * VOLUME * sizeof(Complex));
    initInputFields<<<(VOLUME + 255) / 256, 256>>>(devFermiField, devGaugeField, devStates);

    // Create FermiField and GaugeField objects
    FermiField fermiField(devFermiField, NUM_PARITIES);
    GaugeField gaugeField(devGaugeField, NUM_PARITIES);
    FermiField resultField(devResultField, NUM_PARITIES);

    // Run the stable Conjugate Gradient algorithm
    conjugateGradient(fermiField, gaugeField, resultField, 1000, 1e-6);

    // Copy results from device to host
    Complex *resultFieldHost = new Complex[NUM_PARITIES * VOLUME];
    cudaMemcpy(resultFieldHost, devResultField, NUM_PARITIES * VOLUME * sizeof(Complex), cudaMemcpyDeviceToHost);

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
                        Complex value = resultFieldHost[p * VOLUME + x + DIM_X * (y + DIM_Y * (z + DIM_Z * t))];
                        std::cout << "Result[" << p << "][" << x << "][" << y << "][" << z << "][" << t
                                  << "]: (" << value.real << ", " << value.imag << ")" << std::endl;
                    }
                }
            }
        }
    }

    // Clean up
    cudaFree(devStates);
    cudaFree(devFermiField);
    cudaFree(devGaugeField);
    cudaFree(devResultField);
    delete[] resultFieldHost;

    cudaDeviceSynchronize();
    cudaError_t error = cudaGetLastError();
    if (error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
    }

    return 0;
}
