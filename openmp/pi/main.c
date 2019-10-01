#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 2

static long num_steps = 1000;
double step;

int main()
{
    int i, id;
    double x, pi, sum = 0.0;

    step = 1.0 / (double) num_steps;

    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel private(x, i, id)
    {
        id = omp_get_thread_num();
        for (i = id + 1; i <= num_steps; i = i + NUM_THREADS)
        {
            x = (i - 0.5) * step;
            sum += 4.0 / (1.0 + x*x);
        }
    }
    pi = step * sum;

    printf("Pi is: %.8f\n", pi);
    return 0;
}
