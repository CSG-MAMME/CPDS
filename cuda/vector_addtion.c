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
    cudaGetDeviceCount(&num_devices);

    vec_add_kernel<<<ceil(n, 256)>>>(a, b, c, n);

    return 0;
}
