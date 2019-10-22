#include <stdio.h>

__global__ void hello_world()
{
    printf("Hello World!\n");
}

int main()
    {
        int num_devices;
        cudaGetDeviceCount(&num_devices);

        unsigned int d;
        for (d = 0; d < num_devices; d++)
        {
            cudaSetDevice(d);
            hello_world<<<1,1>>>();
        }

        for (d = 0; d < num_devices; d++)
        {
            cudaSetDevice(d);
            cudaThreadSynchronize();
        } 

        return 0;
}
