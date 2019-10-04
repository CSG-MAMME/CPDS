#include <omp.h>
#include <stdio.h>
#include <stdlib.h>

#define NUM_THREADS 2

#define CRITICAL_MODE 0

static long num_steps = 1000;
double step;

int main()
{
    int i, id;
    double x, pi, sum = 0.0;
    double time_start, time_end;

    time_start = omp_get_wtime();
    step = 1.0 / (double) num_steps;
    omp_set_num_threads(NUM_THREADS);
    // Using the for statement
    #pragma omp parallel for private(x) reduction(+:sum)
    for (i = 1; i <= num_steps; i++)
    {
        x = (i - 0.5) * step;
        sum += 4.0 / (1.0 + x*x);
    }
    pi = step * sum;
    time_end = omp_get_wtime();

    printf("Pi is: \t\t%.8f\n", pi);
    printf("Time elapsed: \t%.8f\n", time_end - time_start);
    return 0;
}
