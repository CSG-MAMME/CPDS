#include <math.h>
#include <float.h>
#include <cuda.h>

__global__ void gpu_Heat (float *h, float *g, float *res, int N)
{
    // CSG: The same Jacobi Method we have used several times
    // CSG: Compute i, j from y and x respectively
    int i = blockIdx.y * blockDim.y + threadIdx.y;
    int j = blockIdx.x * blockDim.x + threadIdx.x;
    // CSG: Sanity check, not to access random memory
    if ((0 < i && i < N - 1) && (0 < j && j < N - 1))
    {
        g[i*N + j]= 0.25 * (h[i*N + (j - 1)] + h[i*N + (j + 1)] +
                            h[(i - 1)*N + j] + h[(i + 1)*N + j]);
    }
    res[i*N + j] = (g[i*N + j] - h[i*N + j])
                    * (g[i*N + j] - h[i*N + j]);
}

// CSG: Basing the code on the class slides (sl 7) on reduction.
__global__ void gpu_Reduction(float *res_in, float *res_out)
{
    extern __shared__ float sdata[];
    unsigned int tid = threadIdx.x;
    // CSG: Optimization from Slide 18
    unsigned int i = blockIdx.x*(blockDim.x*2) + threadIdx.x;
    sdata[tid] = res_in[i] + res_in[i + blockDim.x];
    __syncthreads();

    // CSG: We reduce one full row to one value hence:
    //  - Input Size: np * np
    //  - Output Size: np
    // We use the addressing from Slide 15
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
        {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }

    if (tid == 0) res_out[blockIdx.x] = sdata[0];
}
