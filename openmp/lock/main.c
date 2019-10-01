#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 8

void foo()
{
    omp_lock_t lock;

    omp_init_lock(&lock);
    #pragma omp parallel
    {
        omp_set_lock(&lock);
        printf("Hello world from thread: %d\n", omp_get_thread_num());
        omp_unset_lock(&lock);
    }
}

int main()
{
    omp_set_num_threads(NUM_THREADS);
    foo();
    return 0;
}
