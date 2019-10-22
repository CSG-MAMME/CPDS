### Copy tar File from cpds92001

+ Part2 -> OpenMP
+ Part3 -> MPI
+ Part4 -> CUDA
 
+ `relax_jacobi`:
    - #pragma omp parallel for (Choose the scheduler) -> Slide 50
    - diff must be private
    - reduce for sum

+ `relax_redblack`:
    - how man parallel regions do we need?

+ `relax_gauss`:
    - we must take into account the different dependencies
    - As a consequence, we must introduce certain thread synchronizations (explicit)
    - Either flush (slide 43)
    - Or using the task model (with dependencies)
