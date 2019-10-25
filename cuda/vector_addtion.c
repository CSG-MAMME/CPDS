#include <math.h>
#include <stdio.h>

__global__ void vec_add_kernel(float *a_d, float *b_d, float *c_d, int n)
{
    int i = threadIdx.x + blockDim.x * blockIdx.x;
    c_d[i] = a_d[i] + b_d[i];
}

int main()
{
    int num_devices;
    float A[N], B[N], C[N];
    float *dev_a, *dev_b, *dev_c;

    cudaMalloc(&dev_a, N * sizeof(float));
    cudaMalloc(&dev_b, N * sizeof(float));
    cudaMalloc(&dev_c, N * sizeof(float));

    // Initialize A and B

    cudaMemcpy(dev_a, a, N * sizeof(float), cudaMemcpyHostToDevice);
    cudaMemcpy(dev_b, b, N * sizeof(float), cudaMemcpyHostToDevice);

    dim3 DimGrid(ceil(N/256), 1, 1);
    dim3 DimBlock(256, 1, 1);
    vec_add_kernel<<<DimGrid, DimBlock>>>(dev_a, dev_b, dev_c, N);

    // Copy back result
    cudaMemcpy(c, dev_c, N * sizeof(float), cudaMemcpyDeviceToHost);

    cudaFree(dev_a);
    cudaFree(dev_b);
    cudaFree(dev_c);
    return 0;
}
