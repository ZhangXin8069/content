#include <stdio.h>

__global__ void add(int *a, int *b, int *c) {
    int tid = blockIdx.x;
    if(tid < 1000) {
        c[tid] = a[tid] + b[tid];
    }
}

int main() {
    int *a, *b, *c;
    int size = 1000 * sizeof(int);
    cudaMalloc((void **)&a, size);
    cudaMalloc((void **)&b, size);
    cudaMalloc((void **)&c, size);

    for(int i = 0; i < 1000; i++) {
        a[i] = i;
        b[i] = i * i;
    }

    add<<<1000, 1>>>(a, b, c);

    cudaFree(a);
    cudaFree(b);
    cudaFree(c);

    return 0;
}
