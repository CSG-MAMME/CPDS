==============================================================================
                            CPDS - PAR: OpenMP
==============================================================================

### Loop Parallelism in OpenMP

+ **Worksharing:** divide execution of a code region among members ot a team.
    + Threads _cooperate_ to do some work
    + Two (used) types: `loop` + `single`
    + Lower overhead than `tasks`

+ **Loop Parallelism:** the `for` construct
    ```
    # pragma omp parallel for
    for(init; test; incr)
    ```
    + Iterations must be independent (number of iterations must be computable)
    + Default data-sharing attribute is _shared_
    + The `schedule` clause determines what is executed by each thread:
        + `static[,chunk]`:
        + `dynamic[,chunk]`:
        + `guided[,chunk]`:

