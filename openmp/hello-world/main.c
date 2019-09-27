#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

int main()
{
    //#ifdef _OPENMP
    #pragma omp parallel num_threads(4)
    {
        printf("Hello World from processor: %d\n", omp_get_thread_num());
    }
    return 0;
}
